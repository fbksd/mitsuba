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
#include "mitsuba/render/renderserveradapter.h"
#include "Benchmark/RenderingServer/RenderingServer.h"
#include <QEventLoop>

#define BENCHMARK_SERVER_ON


/**
 *  Discret Probability Density Function sampler.
 *
 *  Permite gerar amostras de acordo com uma pdf por partes dada como entrada no construtor.
 *  Cada valor retornado em sampler() está no intervalo [ 0, fv.size() ).
 *  O custo de sample() é \f$ O(log(n)) \f$ no tamanho de fv e o construtor tem custo linear.
 */
class DPDF
{
    friend class DPDF2;
public:
    DPDF(const float* fv, unsigned int n):
        m_pdf(n),
        m_cpdf(n)
    {
        float sum = std::accumulate(fv, fv + n, 0.f);
        m_funcInt = sum ;

        float nf = 1.f / sum;
        sum = 0.f;

        for(unsigned int i=0; i < n; ++i)
        {
            m_pdf[i] = fv[i] * nf;
            m_cpdf[i] = sum + m_pdf[i];
            sum += m_pdf[i];
        }

        //Garante a última entrada de m_cpdf será 1.f. Pode não ser devido a erro numérico.
        m_cpdf.back() = 1.f;
    }

    DPDF(const std::vector<float>& fv):
        m_pdf(fv.size()),
        m_cpdf(fv.size())
    {
        float sum = accumulate(fv.begin(), fv.end(), 0.f);
        m_funcInt = sum ;

        float nf = 1.f / sum;
        sum = 0.f;

        for(unsigned int i=0; i < fv.size(); ++i)
        {
            m_pdf[i] = fv[i] * nf;
            m_cpdf[i] = sum + m_pdf[i];
            sum += m_pdf[i];
        }

        //Garante a última entrada de m_cpdf será 1.f. Pode não ser devido a erro numérico.
        m_cpdf.back() = 1.f;
    }

    //! Gera um índice aleatório.
    unsigned int sample(float r, float* pdf = NULL)
    {
        float *p = std::lower_bound(&(*m_cpdf.begin()), &(*m_cpdf.end()), r);
        unsigned int i = p - &(*m_cpdf.begin());

        if(pdf && m_pdf.size()) *pdf = m_pdf[i];

        return i;
    }

private:
    std::vector<float> m_pdf;
    std::vector<float> m_cpdf;
    float m_funcInt;
};


class DPDF2
{
public:
    DPDF2(const float* fv, int nx, int ny)
    {
        m_dpdfX.reserve(nx);

        for(int x = 0; x < nx; ++x)
            m_dpdfX.push_back(new DPDF(&fv[x*ny], ny));

        std::vector<float> marginalFunc;
        marginalFunc.reserve(nx);

        for(int x = 0; x < nx; ++x)
            marginalFunc.push_back(m_dpdfX[x]->m_funcInt);

        m_marginal = new DPDF(marginalFunc);
    }

    ~DPDF2()
    {
        for(unsigned int i=0; i< m_dpdfX.size(); ++i)
            delete m_dpdfX[i];

        delete m_marginal;
    }

    //! Gera um índice aleatório.
    void sample(float u, float v, int* x, int* y, float* pdf = NULL)
    {
        float pdfs[2];

        *y = m_marginal->sample(u, &pdfs[1]);
        *x = m_dpdfX[*y]->sample(v, &pdfs[0]);

        if(pdf) *pdf = pdfs[0] * pdfs[1];
    }

private:
    int m_nx, m_ny;
    std::vector<DPDF*> m_dpdfX;
    DPDF* m_marginal;
};



class PixelSampler
{
public:
    PixelSampler(int beginx, int endx, int beginy, int endy, int n, bool isSPP)
    {
        m_beginX = beginx;
        m_endX = endx;
        m_beginY = beginy;
        m_endY = endy;
        m_width = m_endX - m_beginX;
        m_height = m_endY - m_beginY;
        m_xPos = m_beginX;
        m_yPos = m_beginY;
        m_sparseSampleIndex = 0;
        m_random = new Random();

        if(isSPP)
        {
            m_spp = n;
            m_sparseSampleCount = 0;
        }
        else
        {
            m_spp = n / (float)(m_width*m_height);
            m_sparseSampleCount = n % (m_width*m_height);
        }

        if(m_sparseSampleCount)
            m_stage = SPARSE;
        if(m_spp)
            m_stage = SPP;
    }

    virtual bool nextPixel(Point2i* pos)
    {
        int x = m_xPos, y = m_yPos;

        switch(m_stage)
        {
            case SPP:
                if(m_yPos == m_endY || m_xPos == m_endX)
                {
                    m_stage = SPARSE;
                    goto SPARSE_CASE;
                }

                if(++m_xPos == m_endX)
                {
                    m_xPos = m_beginX;
                    ++m_yPos;
                }
                break;
            case SPARSE:
                SPARSE_CASE:
                if(m_sparseSampleIndex++ == m_sparseSampleCount)
                    return false;
                x = m_beginX + m_random->nextUInt(m_width);
                y = m_beginY + m_random->nextUInt(m_height);
        }

        pos->x = x;
        pos->y = y;
        return true;
    }

protected:
    enum Stage
    {
        SPP,
        SPARSE
    };

    int m_beginX;
    int m_endX;
    int m_beginY;
    int m_endY;
    int m_width;
    int m_height;
    int m_xPos, m_yPos;
    int m_spp;
    int m_sparseSampleCount;
    int m_sparseSampleIndex;
    Stage m_stage;
    ref<Random> m_random;
};

class AdaptivePixelSampler: public PixelSampler
{
public:
    AdaptivePixelSampler(int beginx, int endx, int beginy, int endy, int n, bool isSPP, const float* pdf):
        PixelSampler(beginx, endx, beginy, endy, n, isSPP),
        dpdf(pdf, endx - beginx, endy - beginy)
    {}

    virtual bool nextPixel(Point2i* pos)
    {
        bool flag = PixelSampler::nextPixel(pos);
        int x = 0, y = 0;
        dpdf.sample(pos->x, pos->y, &x, &y);
        pos->x = x + m_beginX;
        pos->y = y + m_beginY;
        return flag;
    }

private:
    DPDF2 dpdf;
};


MTS_NAMESPACE_BEGIN

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
    inBuffer = outBuffer = nullptr;
}

SamplingIntegrator::SamplingIntegrator(Stream *stream, InstanceManager *manager)
 : Integrator(stream, manager)
{
    inBuffer = outBuffer = nullptr;
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

#ifdef BENCHMARK_SERVER_ON
    ref<Scheduler> sched = Scheduler::getInstance();
    m_scene =  static_cast<Scene *>(sched->getResource(sceneResID));
    m_sensor = static_cast<Sensor *>(sched->getResource(sensorResID));
    m_originalSampler = static_cast<Sampler *>(sched->getResource(samplerResID, 0));

    QEventLoop eventLoop;
    m_server = new RenderServerAdapter(this);
    QObject::connect(m_server, &RenderServerAdapter::finishRender, &eventLoop, &QEventLoop::quit);
    eventLoop.exec();
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
		const bool &stop, const std::vector< TPoint2<uint8_t> > &points) const {

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

	for (size_t i = 0; i<points.size(); ++i) {
		Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
		if (stop)
			break;

		sampler->generate(offset);

		for (size_t j = 0; j<sampler->getSampleCount(); j++) {
			rRec.newQuery(queryType, sensor->getMedium());
			Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

			if (needsApertureSample)
				apertureSample = rRec.nextSample2D();
			if (needsTimeSample)
				timeSample = rRec.nextSample1D();

			Spectrum spec = sensor->sampleRayDifferential(
				sensorRay, samplePos, apertureSample, timeSample);

			sensorRay.scaleDifferential(diffScaleFactor);

			spec *= Li(sensorRay, rRec);
			block->put(samplePos, spec, rRec.alpha);
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
    info->set<bool>("has_motionblur", shutterTime > 0.0001f);

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

void SamplingIntegrator::evaluateSamples(bool isSPP, int numSamples, int* resultSize)
{
    auto sensorSize = m_sensor->getFilm()->getSize();
    PixelSampler pixelSampler(0, sensorSize.x, 0, sensorSize.y, numSamples, isSPP);

    int totalNumSamples = isSPP ? sensorSize.x * sensorSize.y * numSamples : numSamples;
    *resultSize = totalNumSamples;

    Properties props("MixSampler");
    props.setInteger("sampleCount", numSamples);
    props.setBoolean("isSPP", isSPP);
    props.setInteger("width", sensorSize.x);
    props.setInteger("height", sensorSize.y);
    ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));

    SamplesPipe pipe;
    render(sampler.get(), &pixelSampler, pipe);
}

void SamplingIntegrator::evaluateSamples(bool isSPP, int numSamples, const CropWindow& crop, int* resultSize)
{
    int w = crop.endX - crop.beginX;
    int h = crop.endY - crop.beginY;
    PixelSampler pixelSampler(crop.beginX, crop.endX, crop.beginY, crop.endY, numSamples, isSPP);

    int totalNumSamples = isSPP ? w * h * numSamples : numSamples;
    *resultSize = totalNumSamples;

    Properties props("MixSampler");
    props.setInteger("sampleCount", numSamples);
    props.setBoolean("isSPP", isSPP);
    props.setInteger("width", w);
    props.setInteger("height", h);
    ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));

    SamplesPipe pipe;
    render(sampler.get(), &pixelSampler, pipe);
}

void SamplingIntegrator::evaluateSamples(bool isSPP, int numSamples, const float* pdf, int* resultSize)
{
    auto sensorSize = m_sensor->getFilm()->getSize();
    AdaptivePixelSampler pixelSampler(0, sensorSize.x, 0, sensorSize.y, numSamples, isSPP, pdf);

    int totalNumSamples = isSPP ? sensorSize.x * sensorSize.y * numSamples : numSamples;
    *resultSize = totalNumSamples;

    Properties props("MixSampler");
    props.setInteger("sampleCount", numSamples);
    props.setBoolean("isSPP", isSPP);
    props.setInteger("width", sensorSize.x);
    props.setInteger("height", sensorSize.y);
    ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));

    SamplesPipe pipe;
    render(sampler.get(), &pixelSampler, pipe);
}

void SamplingIntegrator::render(Sampler* sampler, PixelSampler* pixelSampler, SamplesPipe& pipe)
{
    Float diffScaleFactor = 1.0f / std::sqrt((Float) sampler->getSampleCount());

    bool needsApertureSample = m_sensor->needsApertureSample();
    bool needsTimeSample = m_sensor->needsTimeSample();

    RadianceQueryRecord rRec(m_scene, sampler);
    Point2 apertureSample(0.5f);
    Float timeSample = 0.5f;
    RayDifferential sensorRay;

    uint32_t queryType = RadianceQueryRecord::ESensorRay;

    if (!m_sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we don't have to */
        queryType &= ~RadianceQueryRecord::EOpacity;

    Point2i offset;
    while(pixelSampler->nextPixel(&offset))
    {
        sampler->generate(offset);

        for (size_t j = 0; j<sampler->getSampleCount(); j++) {
            SampleBuffer sampleBuffer = pipe.getBuffer();

            rRec.newQuery(queryType, m_sensor->getMedium());
            Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

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

            spec *= Li(sensorRay, rRec);

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
