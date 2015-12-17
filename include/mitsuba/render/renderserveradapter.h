#ifndef RENDERSERVERADAPTER_H
#define RENDERSERVERADAPTER_H

#include <QObject>

class SceneInfo;
class RenderingServer;
class SampleLayout;
class CropWindow;
namespace mitsuba {
    class SamplingIntegrator;
}
using namespace mitsuba;

class RenderServerAdapter: public QObject
{
    Q_OBJECT
public:
    RenderServerAdapter(SamplingIntegrator* integrator);
    ~RenderServerAdapter();

    void setSampleBuffers(float* input, float*output);

signals:
    void finishRender();

private slots:
    void setMaxSPP(int maxSPP);
    void getSceneInfo(SceneInfo* scene);
    void setSampleLayout(const SampleLayout& layout);
    void evaluateSamples(bool isSPP, int numSamples, int* resultSize);
    void evaluateSamplesCrop(bool isSPP, int numSamples, const CropWindow& crop, int* resultSize);
    void evaluateSamplesPDF(bool isSPP, int numSamples, const float* pdf, int* resultSize);

private:
    RenderingServer* server;
    SamplingIntegrator* integrator;
};

#endif // RENDERSERVERADAPTER_H
