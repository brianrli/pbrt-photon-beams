#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_PhotonBeamMAP_H
#define PBRT_INTEGRATORS_PhotonBeamMAP_H

#include "volume.h"
#include "integrator.h"

//Photon Structs
struct PhotonBeam;
struct beamisect;

//BBH Structs
struct BBHBuildNode;

// BBHAccel Forward Declarations
struct BBHPrimitiveInfo;
struct LinearBBHNode;
class BBHAccel;

// PhotonIntegrator Declarations
class PhotonBeamIntegrator : public VolumeIntegrator {

public:
    // PhotonIntegrator Public Methods
    PhotonBeamIntegrator(int nbeams, float ssize);
    ~PhotonBeamIntegrator(){}
    
    Spectrum Li(const Scene *scene,
                const Renderer *renderer,
                const RayDifferential &ray,
                const Sample *sample,
                RNG &rng,
                Spectrum *transmittance,
                MemoryArena &arena)const;
    
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
    float stepSize;
    float blur; 
    
    // Declare sample parameters for light source sampling
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;
    BSDFSampleOffsets bsdfGatherSampleOffsets, indirGatherSampleOffsets;
    int nVolumePaths;
    
    //BBH Hierarchy
    BBHAccel *BeamMap;
};

class PhotonBeamShootingTask : public Task {
public:
    PhotonBeamShootingTask(int tn, float ti, Mutex &m, PhotonBeamIntegrator *in,
                           ProgressReporter &prog, bool &at, uint32_t &ns,
                           vector<PhotonBeam> &vol, Distribution1D *distrib, const Scene *sc,
                           const Renderer *sr)
    : taskNum(tn), time(ti), mutex(m), integrator(in), progress(prog),
    abortTasks(at), nshot(ns), lightDistribution(distrib), scene(sc), renderer (sr),
    PhotonBeams(vol){}
    
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

class BBHAccel
{
public:
    // VBVHAccel Public Methods
    BBHAccel(const vector<PhotonBeam> &p, uint32_t mp);
    ~BBHAccel();
    
    //Intersect
    bool Intersect(const Ray &ray, vector<beamisect> &intersections) const;
    
    BBox WorldBound() const;
    //bool CanIntersect() const { return true; }
private:
    // VBVHAccel Private Methods
    BBHBuildNode *recursiveBuild(MemoryArena &buildArena,
                                 vector<BBHPrimitiveInfo> &buildData, uint32_t start, uint32_t end,
                                 uint32_t *totalNodes, vector<PhotonBeam> &orderedBeams);
    uint32_t flattenBBHTree(BBHBuildNode *node, uint32_t *offset);
    
    // VBVHAccel Private Data
    uint32_t maxBeamsInNode;
    vector<PhotonBeam> beams;
    LinearBBHNode *nodes;
};

BBHAccel *CreateBBHAccelerator(const vector<PhotonBeam> &beams);

PhotonBeamIntegrator *CreatePhotonBeamIntegrator(const ParamSet &params);


#endif // PBRT_INTEGRATORS_PHOTONMAP_H
