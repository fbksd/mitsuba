
#include <mitsuba/render/sampler.h>
#include <mitsuba/core/plugin.h>

MTS_NAMESPACE_BEGIN


class MixSampler : public Sampler {
public:
    MixSampler() : Sampler(Properties()) { }

    MixSampler(const Properties &props) : Sampler(props) {
        /* Number of samples per pixel when used with a sampling-based integrator */
        int n = props.getSize("sampleCount", 4);
        m_mainSampler = m_sparseSampler = nullptr;
        bool isSPP = props.getBoolean("isSPP", true);
        int w = props.getInteger("width", 1);
        int h = props.getInteger("height", 1);
        m_numPixels = w * h;
        m_pixelIndex = 0;
        m_sparseSampleCount = 0;

        if(isSPP)
            m_spp = n;
        else
        {
            m_spp = n / (float)(w*h);
            m_sparseSampleCount = n % (w*h);
        }

        if(m_sparseSampleCount)
        {
            Properties p("independent");
            p.setInteger("sampleCount", 1);
            m_sparseSampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), p));
            m_stage = SPARSE;
            m_sampleCount = 1;
        }
        if(m_spp)
        {
            Properties p("independent");
            p.setInteger("sampleCount", m_spp);
            m_mainSampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), p));
            m_stage = MAIN;
            m_sampleCount = m_spp;
        }
    }

    MixSampler(Stream *stream, InstanceManager *manager)
     : Sampler(stream, manager) {
        m_mainSampler = static_cast<Sampler *>(manager->getInstance(stream));
        m_sparseSampler = static_cast<Sampler *>(manager->getInstance(stream));
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Sampler::serialize(stream, manager);
        manager->serialize(stream, m_mainSampler.get());
        manager->serialize(stream, m_sparseSampler.get());
    }

    ref<Sampler> clone() {
        ref<MixSampler> sampler = new MixSampler();
        sampler->m_sampleCount = m_sampleCount;
        sampler->m_mainSampler = m_mainSampler->clone();
        sampler->m_sparseSampler = m_sparseSampler->clone();
        return sampler.get();
    }

    void generate(const Point2i &p) {
        switch(m_stage)
        {
            case MAIN:
                if(m_pixelIndex++ < m_numPixels)
                {
                    m_mainSampler->generate(p);
                    break;
                }
                else
                {
                    m_stage = SPARSE;
                    m_sampleCount = 1;
                }

            case SPARSE:
                m_sparseSampler->generate(p);
        }
    }

    Float next1D() {
        switch(m_stage)
        {
            case MAIN:
                return m_mainSampler->next1D();
            case SPARSE:
                return m_sparseSampler->next1D();
        }
    }

    Point2 next2D() {
        switch(m_stage)
        {
            case MAIN:
                return m_mainSampler->next2D();
            case SPARSE:
                return m_sparseSampler->next2D();
        }
    }

    void advance()
    {
        Sampler::advance();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "MixSampler[" << endl
            << "  sampleCount = " << m_sampleCount << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    enum Stage
    {
        MAIN,
        SPARSE
    };

    int m_spp;
    int m_sparseSampleCount;
    int m_numPixels;
    int m_pixelIndex;
    Stage m_stage;
    ref<Sampler> m_mainSampler;
    ref<Sampler> m_sparseSampler;
};


MTS_IMPLEMENT_CLASS_S(MixSampler, false, Sampler)
MTS_EXPORT_PLUGIN(MixSampler, "Mix Sampler");
MTS_NAMESPACE_END
