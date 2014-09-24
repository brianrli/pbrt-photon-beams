#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_PhotonBeamMAP_H
#define PBRT_INTEGRATORS_PhotonBeamMAP_H

#include "volume.h"
#include "integrator.h"

//Photon Structs
struct PhotonBeam;

// PhotonIntegrator Declarations
class PhotonBeamIntegrator : public VolumeIntegrator {

public:
    // PhotonIntegrator Public Methods
    
    PhotonBeamIntegrator(int nbeams);
    
    Spectrum Li(const Scene *scene,
                const Renderer *renderer,
                const RayDifferential &ray,
                const Sample *sample,
                RNG &rng,
                Spectrum *transmittance,
                MemoryArena &arena) const;
    
    Spectrum Transmittance(const Scene *scene,
                           const Renderer *,
                           const RayDifferential &ray,
                           const Sample *sample,
                           RNG &rng,
                           MemoryArena &arena) const;
    
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    
    void Preprocess(const Scene *scene, const Camera *camera, const Renderer *renderer);
    
private:
    // PhotonIntegrator Private Methods
    friend class PhotonBeamShootingTask;
    
    // PhotonIntegrator Private Data
    uint32_t nPhotonBeamsWanted;
    int tauSampleOffset, scatterSampleOffset;
    float blur; 
    
    // Declare sample parameters for light source sampling
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;
    BSDFSampleOffsets bsdfGatherSampleOffsets, indirGatherSampleOffsets;
    int nVolumePaths;
};

class PhotonBeamShootingTask : public Task {
public:
    PhotonBeamShootingTask(int tn, float ti, Mutex &m, PhotonBeamIntegrator *in,
                           ProgressReporter &prog, bool &at, uint32_t &ns,
                           vector<PhotonBeam> &vol, Distribution1D *distrib, const Scene *sc,
                           const Renderer *sr)
    : taskNum(tn), time(ti), mutex(m), integrator(in), progress(prog),
    abortTasks(at), nshot(ns), lightDistribution(distrib), scene(sc), renderer (sr),
    PhotonBeams(vol){ }
    
    void Run();
    
    int taskNum;
    float time;
    Mutex &mutex;
    PhotonBeamIntegrator *integrator;
    ProgressReporter &progress;
    bool &abortTasks;
    vector<PhotonBeam> &PhotonBeams;
    uint32_t &nshot;
    const Distribution1D *lightDistribution;
    const Scene *scene;
    const Renderer *renderer;
};



PhotonBeamIntegrator *CreatePhotonBeamIntegrator(const ParamSet &params);


#endif // PBRT_INTEGRATORS_PHOTONMAP_H
