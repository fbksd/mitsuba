#include "mitsuba/render/renderserveradapter.h"
#include "Benchmark/RenderingServer/RenderingServer.h"
#include "mitsuba/render/renderjob.h"

RenderServerAdapter::RenderServerAdapter(SamplingIntegrator* integrator):
    QObject(NULL),
    server(new RenderingServer),
    integrator(integrator)
{
    QObject::connect(server, &RenderingServer::getSceneInfo, this, &RenderServerAdapter::getSceneInfo);
    QObject::connect(server, &RenderingServer::evaluateSamples, this, &RenderServerAdapter::evaluateSamples);
    QObject::connect(server, &RenderingServer::evaluateSamplesCrop, this, &RenderServerAdapter::evaluateSamplesCrop);
    QObject::connect(server, &RenderingServer::evaluateSamplesPDF, this, &RenderServerAdapter::evaluateSamplesPDF);
    QObject::connect(server, &RenderingServer::finishRender, this, &RenderServerAdapter::finishRender);
    server->startServer(2227);
}

void RenderServerAdapter::getSceneInfo(SceneInfo *scene)
{
    integrator->getSceneInfo(scene);
}

void RenderServerAdapter::evaluateSamples(bool isSPP, int numSamples, int *resultSize)
{
    integrator->evaluateSamples(isSPP, numSamples, resultSize);
}

void RenderServerAdapter::evaluateSamplesCrop(bool isSPP, int numSamples, const CropWindow &crop, int *resultSize)
{
    integrator->evaluateSamples(isSPP, numSamples, crop, resultSize);
}

void RenderServerAdapter::evaluateSamplesPDF(bool isSPP, int numSamples, const float *pdf, int *resultSize)
{
    integrator->evaluateSamples(isSPP, numSamples, pdf, resultSize);
}

