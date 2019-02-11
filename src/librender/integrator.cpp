/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/statistics.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/renderproc.h>

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/random.h>
#include <fbksd/renderer/RenderingServer.h>
#include <iostream>

#define BENCHMARK_SERVER_ON

MTS_NAMESPACE_BEGIN


class PixelSampler
{
public:
    PixelSampler(int beginx, int endx, int beginy, int endy, int n)
    {
        m_beginX = beginx;
        m_endX = endx;
        m_beginY = beginy;
        m_endY = endy;
        m_width = m_endX - m_beginX;
        m_height = m_endY - m_beginY;
        m_sampleIndex = 0;
        m_random = new Random();
        m_numSamples = n;
    }

    virtual bool nextPixel(Point2i* pos)
    {
        if(m_sampleIndex++ == m_numSamples)
            return false;
        pos->x = m_beginX + m_random->nextUInt(m_width);
        pos->y = m_beginY + m_random->nextUInt(m_height);
        return true;
    }

protected:
    int m_beginX;
    int m_endX;
    int m_beginY;
    int m_endY;
    int m_width;
    int m_height;
    int m_numSamples;
    int m_sampleIndex;
    ref<Random> m_random;
};


Integrator::Integrator(const Properties &props)
 : NetworkedObject(props) { }

Integrator::Integrator(Stream *stream, InstanceManager *manager)
 : NetworkedObject(stream, manager) { }

bool Integrator::preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int sensorResID, int samplerResID) { return true; }
void Integrator::postprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int sensorResID, int samplerResID) { }
void Integrator::serialize(Stream *stream, InstanceManager *manager) const {
	NetworkedObject::serialize(stream, manager);
}
void Integrator::configureSampler(const Scene *scene, Sampler *sampler) {
	/* Prepare the sampler for bucket-based rendering */
	sampler->setFilmResolution(scene->getFilm()->getCropSize(),
		getClass()->derivesFrom(MTS_CLASS(SamplingIntegrator)));
}
const Integrator *Integrator::getSubIntegrator(int idx) const { return NULL; }

SamplingIntegrator::SamplingIntegrator(const Properties &props)
 : Integrator(props)
{
}

SamplingIntegrator::SamplingIntegrator(Stream *stream, InstanceManager *manager)
 : Integrator(stream, manager)
{
}

void SamplingIntegrator::serialize(Stream *stream, InstanceManager *manager) const {
	Integrator::serialize(stream, manager);
}

Spectrum SamplingIntegrator::E(const Scene *scene, const Intersection &its,
		const Medium *medium, Sampler *sampler, int nSamples, bool handleIndirect) const {
	Spectrum E(0.0f);
	RadianceQueryRecord query(scene, sampler);
	DirectSamplingRecord dRec(its);
	Frame frame(its.shFrame.n);

	sampler->generate(Point2i(0));
	for (int i=0; i<nSamples; i++) {
		/* Sample the direct illumination component */
		int maxIntermediateInteractions = -1;
		Spectrum directRadiance = scene->sampleAttenuatedEmitterDirect(
			dRec, its, medium, maxIntermediateInteractions, query.nextSample2D());

		if (!directRadiance.isZero()) {
			Float dp = dot(dRec.d, its.shFrame.n);
			if (dp > 0)
				E += directRadiance * dp;
		}

		/* Sample the indirect illumination component */
		if (handleIndirect) {
			query.newQuery(RadianceQueryRecord::ERadianceNoEmission, medium);
			Vector d = frame.toWorld(warp::squareToCosineHemisphere(query.nextSample2D()));
			++query.depth;
			query.medium = medium;
			E += Li(RayDifferential(its.p, d, its.time), query) * M_PI;
		}

		sampler->advance();
	}

	return E / (Float) nSamples;
}

void SamplingIntegrator::cancel() {
	if (m_process)
		Scheduler::getInstance()->cancel(m_process);
}

bool SamplingIntegrator::render(Scene *scene,
		RenderQueue *queue, const RenderJob *job,
		int sceneResID, int sensorResID, int samplerResID) {

    m_queue = queue;
    m_job = const_cast<RenderJob*>(job);
    m_sceneResID = sceneResID;
    m_sensorResID = sensorResID;
    m_samplerResID = samplerResID;

#ifdef BENCHMARK_SERVER_ON
    ref<Scheduler> sched = Scheduler::getInstance();
    m_scene =  static_cast<Scene *>(sched->getResource(sceneResID));
    m_sensor = static_cast<Sensor *>(sched->getResource(sensorResID));
    m_originalSampler = static_cast<Sampler *>(sched->getResource(samplerResID, 0));

    RenderingServer server;
    server.onGetTileSize([](){ return 32; });
    server.onSetParameters([this](const SampleLayout& layout){
        m_layout = layout;
    });
    server.onGetSceneInfo([this](){
        SceneInfo info;
        this->getSceneInfo(&info);
        return info;
    });
    server.onEvaluateSamples([this](int64_t spp, int64_t remainingCount, int pipeSize){
        this->evaluateSamples(spp, remainingCount, pipeSize);
    });
    server.onLastTileConsumed([this](){
        if(m_process)
        {
            Scheduler::getInstance()->wait(m_process);
            m_process = nullptr;
        }
        if(m_sparseProcess)
        {
            Scheduler::getInstance()->wait(m_sparseProcess);
            m_sparseProcess = nullptr;
        }
        Scheduler::getInstance()->unregisterResource(m_integratorResID);
    });
    server.run();
    return true;
#else
    ref<Scheduler> sched = Scheduler::getInstance();
    ref<Sensor> sensor = static_cast<Sensor *>(sched->getResource(sensorResID));
    ref<Film> film = sensor->getFilm();
    size_t nCores = sched->getCoreCount();
    const Sampler *sampler = static_cast<const Sampler *>(sched->getResource(samplerResID, 0));
    size_t sampleCount = sampler->getSampleCount();

    Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SIZE_T_FMT
        " %s, " SSE_STR ") ..", film->getCropSize().x, film->getCropSize().y,
        sampleCount, sampleCount == 1 ? "sample" : "samples", nCores,
        nCores == 1 ? "core" : "cores");

    /* This is a sampling-based integrator - parallelize */
    ref<ParallelProcess> proc = new BlockedRenderProcess(job,
        queue, scene->getBlockSize());
    int integratorResID = sched->registerResource(this);
    proc->bindResource("integrator", integratorResID);
    proc->bindResource("scene", sceneResID);
    proc->bindResource("sensor", sensorResID);
    proc->bindResource("sampler", samplerResID);
    scene->bindUsedResources(proc);
    bindUsedResources(proc);
    sched->schedule(proc);

    m_process = proc;
    sched->wait(proc);
    m_process = NULL;
    sched->unregisterResource(integratorResID);

    return proc->getReturnStatus() == ParallelProcess::ESuccess;
#endif
}

void SamplingIntegrator::bindUsedResources(ParallelProcess *) const {
	/* Do nothing by default */
}

void SamplingIntegrator::wakeup(ConfigurableObject *parent,
	std::map<std::string, SerializableObject *> &) {
	/* Do nothing by default */
}

void SamplingIntegrator::renderBlock(const Scene *scene,
		const Sensor *sensor, Sampler *sampler, ImageBlock *block,
        const bool &stop, const std::vector< TPoint2<uint8_t> > &points, bool seekPipeByPixel) const {

	Float diffScaleFactor = 1.0f /
		std::sqrt((Float) sampler->getSampleCount());

	bool needsApertureSample = sensor->needsApertureSample();
	bool needsTimeSample = sensor->needsTimeSample();

	RadianceQueryRecord rRec(scene, sampler);
	Point2 apertureSample(0.5f);
	Float timeSample = 0.5f;
	RayDifferential sensorRay;

	block->clear();

	uint32_t queryType = RadianceQueryRecord::ESensorRay;

	if (!sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we don't have to */
		queryType &= ~RadianceQueryRecord::EOpacity;

#ifdef BENCHMARK_SERVER_ON
    const auto& blockOffset = Vector2i(block->getOffset());
    SamplesPipe pipe({blockOffset.x, blockOffset.y},
                     {blockOffset.x + block->getWidth(), blockOffset.y + block->getHeight()},
                     points.size() * sampler->getSampleCount());
#endif

	for (size_t i = 0; i<points.size(); ++i) {
		Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
        if (stop)
			break;

		sampler->generate(offset);

#ifdef BENCHMARK_SERVER_ON
        if(seekPipeByPixel)
            pipe.seek(offset.x, offset.y);
#endif

		for (size_t j = 0; j<sampler->getSampleCount(); j++) {
#ifdef BENCHMARK_SERVER_ON
            SampleBuffer sampleBuffer = pipe.getBuffer();
#endif
			rRec.newQuery(queryType, sensor->getMedium());
			Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

			if (needsApertureSample)
				apertureSample = rRec.nextSample2D();
			if (needsTimeSample)
				timeSample = rRec.nextSample1D();

#ifdef BENCHMARK_SERVER_ON
            // Sometimes sample spill out of the correct pixel due to rounding errors.
            // Trying to fix this here.
            int x = samplePos.x;
            int y = samplePos.y;
            if(x > offset.x)
                samplePos.x = (offset.x + 1) - Epsilon;
            if(y > offset.y)
                samplePos.y = (offset.y + 1) - Epsilon;

            samplePos.x = sampleBuffer.set(IMAGE_X, samplePos.x);
            samplePos.y = sampleBuffer.set(IMAGE_Y, samplePos.y);
            apertureSample.x = sampleBuffer.set(LENS_U, apertureSample.x);
            apertureSample.y = sampleBuffer.set(LENS_V, apertureSample.y);
            timeSample = sampleBuffer.set(TIME, timeSample);
#endif

            Spectrum spec = sensor->sampleRayDifferential(
				sensorRay, samplePos, apertureSample, timeSample);

			sensorRay.scaleDifferential(diffScaleFactor);

#ifdef BENCHMARK_SERVER_ON
            spec *= Li(sensorRay, rRec, &sampleBuffer);
#else
            spec *= Li(sensorRay, rRec, nullptr);
#endif
            bool validSample = block->put(samplePos, spec, rRec.alpha);
			sampler->advance();

#ifdef BENCHMARK_SERVER_ON
            float r = 0.f, g = 0.f, b = 0.f;
            if(validSample)
                spec.toLinearRGB(r, g, b);
            sampleBuffer.set(COLOR_R, r);
            sampleBuffer.set(COLOR_G, g);
            sampleBuffer.set(COLOR_B, b);
            pipe << sampleBuffer;
#endif
		}
    }
}

void SamplingIntegrator::renderSamples(const Scene *scene, const Sensor *sensor, Sampler *sampler, int numSamples) const
{
    const auto sensorSize = sensor->getFilm()->getSize();
    SamplesPipe pipe({0, 0}, {sensorSize.x, sensorSize.y}, numSamples);

    Float diffScaleFactor = 1.0f / std::sqrt((Float) sampler->getSampleCount());
    bool needsApertureSample = m_sensor->needsApertureSample();
    bool needsTimeSample = m_sensor->needsTimeSample();

    RadianceQueryRecord rRec(scene, sampler);
    Point2 apertureSample(0.5f);
    Float timeSample = 0.5f;
    RayDifferential sensorRay;

    uint32_t queryType = RadianceQueryRecord::ESensorRay;

    if (!m_sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we don't have to */
        queryType &= ~RadianceQueryRecord::EOpacity;

    PixelSampler pixelSampler(0, sensorSize.x, 0, sensorSize.y, numSamples);
    Point2i offset;
    while(pixelSampler.nextPixel(&offset))
    {
        sampler->generate(offset);

        for (size_t j = 0; j<sampler->getSampleCount(); j++) {
            SampleBuffer sampleBuffer = pipe.getBuffer();

            rRec.newQuery(queryType, m_sensor->getMedium());
            Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

            // Sometimes sample spill out of the correct pixel due to rounding errors.
            // Trying to fix this here.
            int x = samplePos.x;
            int y = samplePos.y;
            if(x > offset.x)
                samplePos.x = (offset.x + 1) - Epsilon;
            if(y > offset.y)
                samplePos.y = (offset.y + 1) - Epsilon;

            if (needsApertureSample)
                apertureSample = rRec.nextSample2D();
            if (needsTimeSample)
                timeSample = rRec.nextSample1D();

            samplePos.x = sampleBuffer.set(IMAGE_X, samplePos.x);
            samplePos.y = sampleBuffer.set(IMAGE_Y, samplePos.y);
            apertureSample.x = sampleBuffer.set(LENS_U, apertureSample.x);
            apertureSample.y = sampleBuffer.set(LENS_V, apertureSample.y);
            timeSample = sampleBuffer.set(TIME, timeSample);

            Spectrum spec = m_sensor->sampleRayDifferential(
                sensorRay, samplePos, apertureSample, timeSample);

            sensorRay.scaleDifferential(diffScaleFactor);

            spec *= Li(sensorRay, rRec, &sampleBuffer);

            float r, g, b;
            spec.toLinearRGB(r, g, b);
            sampleBuffer.set(COLOR_R, r);
            sampleBuffer.set(COLOR_G, g);
            sampleBuffer.set(COLOR_B, b);
            pipe << sampleBuffer;

            sampler->advance();
        }
    }
}

void SamplingIntegrator::getSceneInfo(SceneInfo *info)
{
    Vector2i filmSize = m_sensor->getFilm()->getSize();
    info->set<int>("width", filmSize.x);
    info->set<int>("height", filmSize.y);
    float shutterOpen = m_sensor->getShutterOpen();
    float shutterTime = m_sensor->getShutterOpenTime();
    info->set<float>("shutter_open", shutterOpen);
    info->set<float>("shutter_close", shutterOpen + shutterTime);
    info->set<bool>("has_motion_blur", shutterTime > 0.0001f);
    size_t sampleCount = m_originalSampler->getSampleCount();
    info->set<int>("max_spp", sampleCount);
    info->set<int>("max_samples", sampleCount*filmSize.x*filmSize.y);

    // TODO:
//    auto shapes = m_scene->getShapes();
//    bool hasAreaLights = false;
//    for(Shape* shape : shapes)
//    {
//        shape->isEmitter()
//        if(light->IsDeltaLight() == false)
//        {
//            hasAreaLights = true;
//            break;
//        }
//    }
//    info->set<bool>("has_area_lights", hasAreaLights);
}

void SamplingIntegrator::evaluateSamples(int64_t spp, int64_t remaining, int pipeSize)
{
    bool hasInputParams = m_layout.hasInput("IMAGE_X") || m_layout.hasInput("IMAGE_Y");

    if(spp)
    {
        ref<Scheduler> sched = Scheduler::getInstance();
        // replace the sampler that comes with the scene by one with spp samples
        sched->unregisterResource(m_samplerResID);
        Properties p;
        if(math::isPowerOfTwo(spp))
            p = Properties("ldsampler");
        else
            p = Properties("independent");
        p.setInteger("sampleCount", spp);
        ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), p));
        std::vector<SerializableObject *> samplers(sched->getCoreCount());
        for (size_t i=0; i<sched->getCoreCount(); ++i) {
            ref<Sampler> clonedSampler = sampler->clone();
            clonedSampler->incRef();
            samplers[i] = clonedSampler.get();
        }
        m_samplerResID = sched->registerMultiResource(samplers);
        for (size_t i=0; i<sched->getCoreCount(); ++i)
            samplers[i]->decRef();
        ref<Sensor> sensor = static_cast<Sensor *>(sched->getResource(m_sensorResID));
        ref<Film> film = sensor->getFilm();
        size_t nCores = sched->getCoreCount();
        size_t sampleCount = sampler->getSampleCount();

        Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SIZE_T_FMT
            " %s, " SSE_STR ") ..", film->getCropSize().x, film->getCropSize().y,
            sampleCount, sampleCount == 1 ? "spp" : "spps", nCores,
            nCores == 1 ? "core" : "cores");

        /* This is a sampling-based integrator - parallelize */
        ref<ParallelProcess> proc = new BlockedRenderProcess(m_job, m_queue, m_scene->getBlockSize(), spp, m_layout.getSampleSize(), !hasInputParams);
        m_integratorResID = sched->registerResource(this);
        proc->bindResource("integrator", m_integratorResID);
        proc->bindResource("scene", m_sceneResID);
        proc->bindResource("sensor", m_sensorResID);
        proc->bindResource("sampler", m_samplerResID);
        m_scene->bindUsedResources(proc);
        bindUsedResources(proc);
        sched->schedule(proc);
        m_process = proc;
    }

    if(remaining)
    {
        ref<Scheduler> sched = Scheduler::getInstance();
        // replace the sampler that comes with the scene by one with spp samples
        sched->unregisterResource(m_samplerResID);
        Properties p("independent");
        p.setInteger("sampleCount", 1);
        ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), p));
        std::vector<SerializableObject *> samplers(sched->getCoreCount());
        for (size_t i=0; i<sched->getCoreCount(); ++i) {
            ref<Sampler> clonedSampler = sampler->clone();
            clonedSampler->incRef();
            samplers[i] = clonedSampler.get();
        }
        m_samplerResID = sched->registerMultiResource(samplers);
        for (size_t i=0; i<sched->getCoreCount(); ++i)
            samplers[i]->decRef();
        ref<Sensor> sensor = static_cast<Sensor *>(sched->getResource(m_sensorResID));
        ref<Film> film = sensor->getFilm();
        size_t nCores = sched->getCoreCount();

        Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SIZE_T_FMT
            " %s, " SSE_STR ") ..", film->getCropSize().x, film->getCropSize().y,
            remaining, "samples", nCores,
            nCores == 1 ? "core" : "cores");

        /* This is a sampling-based integrator - parallelize */
        ref<ParallelProcess> proc = new SparseRenderProcess(m_job, m_queue, remaining, pipeSize);
        if(!spp)
            m_integratorResID = sched->registerResource(this);
        proc->bindResource("integrator", m_integratorResID);
        proc->bindResource("scene", m_sceneResID);
        proc->bindResource("sensor", m_sensorResID);
        proc->bindResource("sampler", m_samplerResID);
        m_scene->bindUsedResources(proc);
        bindUsedResources(proc);
        sched->schedule(proc);
        m_sparseProcess = proc;
    }
}


MonteCarloIntegrator::MonteCarloIntegrator(const Properties &props) : SamplingIntegrator(props) {
	/* Depth to begin using russian roulette */
	m_rrDepth = props.getInteger("rrDepth", 5);

	/* Longest visualized path depth (\c -1 = infinite).
	   A value of \c 1 will visualize only directly visible light sources.
	   \c 2 will lead to single-bounce (direct-only) illumination, and so on. */
	m_maxDepth = props.getInteger("maxDepth", -1);

	/**
	 * This parameter specifies the action to be taken when the geometric
	 * and shading normals of a surface don't agree on whether a ray is on
	 * the front or back-side of a surface.
	 *
	 * When \c strictNormals is set to \c false, the shading normal has
	 * precedence, and rendering proceeds normally at the risk of
	 * introducing small light leaks (this is the default).
	 *
	 * When \c strictNormals is set to \c true, the random walk is
	 * terminated when encountering such a situation. This may
	 * lead to silhouette darkening on badly tesselated meshes.
	 */
	m_strictNormals = props.getBoolean("strictNormals", false);

	/**
	 * When this flag is set to true, contributions from directly
	 * visible emitters will not be included in the rendered image
	 */
	m_hideEmitters = props.getBoolean("hideEmitters", false);

	if (m_rrDepth <= 0)
		Log(EError, "'rrDepth' must be set to a value greater than zero!");

	if (m_maxDepth <= 0 && m_maxDepth != -1)
		Log(EError, "'maxDepth' must be set to -1 (infinite) or a value greater than zero!");
}

MonteCarloIntegrator::MonteCarloIntegrator(Stream *stream, InstanceManager *manager)
	: SamplingIntegrator(stream, manager) {
	m_rrDepth = stream->readInt();
	m_maxDepth = stream->readInt();
	m_strictNormals = stream->readBool();
	m_hideEmitters = stream->readBool();
}

void MonteCarloIntegrator::serialize(Stream *stream, InstanceManager *manager) const {
	SamplingIntegrator::serialize(stream, manager);
	stream->writeInt(m_rrDepth);
	stream->writeInt(m_maxDepth);
	stream->writeBool(m_strictNormals);
	stream->writeBool(m_hideEmitters);
}

std::string RadianceQueryRecord::toString() const {
	std::ostringstream oss;
	oss << "RadianceQueryRecord[" << endl
		<< "  type = { ";
	if (type & EEmittedRadiance) oss << "emitted ";
	if (type & ESubsurfaceRadiance) oss << "subsurface ";
	if (type & EDirectSurfaceRadiance) oss << "direct ";
	if (type & EIndirectSurfaceRadiance) oss << "indirect ";
	if (type & ECausticRadiance) oss << "caustic ";
	if (type & EDirectMediumRadiance) oss << "inscatteredDirect ";
	if (type & EIndirectMediumRadiance) oss << "inscatteredIndirect ";
	if (type & EDistance) oss << "distance ";
	if (type & EOpacity) oss << "opacity ";
	if (type & EIntersection) oss << "intersection ";
	oss << "}," << endl
		<< "  depth = " << depth << "," << endl
		<< "  its = " << indent(its.toString()) << endl
		<< "  alpha = " << alpha << "," << endl
		<< "  extra = " << extra << "," << endl
		<< "]" << endl;
	return oss.str();
}


MTS_IMPLEMENT_CLASS(Integrator, true, NetworkedObject)
MTS_IMPLEMENT_CLASS(SamplingIntegrator, true, Integrator)
MTS_IMPLEMENT_CLASS(MonteCarloIntegrator, true, SamplingIntegrator)
MTS_NAMESPACE_END
