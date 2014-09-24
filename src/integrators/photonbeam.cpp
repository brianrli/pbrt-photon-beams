// integrators/photonmap.cpp*
#include "stdafx.h"
#include "integrators/photonbeam.h"
#include "scene.h"
#include "montecarlo.h"
#include "sampler.h"
#include "progressreporter.h"
#include "intersection.h"
#include "paramset.h"
#include "camera.h"

struct PhotonBeam {
    //constructor
    PhotonBeam(const Point &start,
               const Point &end,
               const Spectrum &wt,//alpha
               const Vector &w, //direction
               float radius)
    : start(start), alpha(wt), dir(w), ri(radius),zmin(0)
    {
        //necessary transformations
        WorldToObject = LookAt(start, end, Cross(end-start,Vector(1,0,0)));
        ObjectToWorld = Inverse(WorldToObject);
        //zmax = ???
    }

    //computes intersect (bastardized cylinder code)
    bool Intersect(Ray &r,float &t0, float &t1){
        float phi;
        Point phit;
        // Transform _Ray_ to object space
        Ray ray;
        WorldToObject(r, &ray);
        
        // Compute quadratic cylinder coefficients
        float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y;
        float B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y);
        float C = ray.o.x*ray.o.x + ray.o.y*ray.o.y - ri*ri;
        
        // Solve quadratic equation for _t_ values
        if (!Quadratic(A, B, C, &t0, &t1))
            return false;
        
        // Compute intersection distance along ray
        if (t0 > ray.maxt || t1 < ray.mint)
            return false;
        float thit = t0;
        if (t0 < ray.mint) {
            thit = t1;
            if (thit > ray.maxt) return false;
        }
        
        // Compute cylinder hit point and $\phi$
        phit = ray(thit);
        phi = atan2f(phit.y, phit.x);
        if (phi < 0.) phi += 2.f*M_PI;
        
        // Test cylinder intersection against clipping parameters
        if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
            if (thit == t1) return false;
            thit = t1;
            if (t1 > ray.maxt) return false;
            // Compute cylinder hit point and $\phi$
            phit = ray(thit);
            phi = atan2f(phit.y, phit.x);
            if (phi < 0.) phi += 2.f*M_PI;
            if (phit.z < zmin || phit.z > zmax || phi > phiMax)
                return false;
        }
        return true;
    }
    
    //returns photon beam at split point. truncates this photon beam
    PhotonBeam split(float split_point){}

    /** Members **/
    Point start,end;
    Spectrum alpha;
    Vector dir;
    Transform ObjectToWorld, WorldToObject;
    float ri;
    BBox bound;
    
    /** Cylinder **/
    float phiMax =2.f*M_PI;
    float zmin, zmax;
};

//get cross product of vector with x axis
//use that as up vector, and bastardize look at function
//get transform, bastardize cylinder transform


//write split function

//BVH implementation

inline bool unsuccessful(uint32_t needed, uint32_t found, uint32_t shot) {
    return (found < needed && (found == 0 || found < shot / 1024));
}

void PhotonBeamShootingTask::Run() {
    
    // Declare local variables for _PhotonBeamShootingTask_
    VolumeRegion *volume = scene->volumeRegion;
    
    MemoryArena arena;
    RNG rng(31 * taskNum);
    
    vector<PhotonBeam> localPhotonBeams; //new data structure
    uint32_t totalPaths = 0;
    PermutedHalton halton(6, rng);

    //set finish conditions
    bool volumeDone = (integrator->nPhotonBeamsWanted == 0);

    while (true) {
        const uint32_t blockSize = 500;
        for (uint32_t i = 0; i < blockSize; ++i) {
            
            float u[6];
            halton.Sample(++totalPaths, u); //random sample
            
            // Choose light to shoot photon from
            float lightPdf;
            int lightNum = lightDistribution->SampleDiscrete(u[0], &lightPdf);
            const Light *light = scene->lights[lightNum]; //choose light
            
            // Generate _photonRay_ from light source and initialize _alpha_
            RayDifferential photonRay;
            float pdf;
            LightSample ls(u[1], u[2], u[3]);
            Normal Nl;

            //Alpha (PBRT)
            Spectrum Le = light->Sample_L(scene, ls, u[4], u[5],
                                          time, &photonRay, &Nl, &pdf);
            if (pdf == 0.f || Le.IsBlack()) continue;
            Spectrum alpha = (AbsDot(Nl, photonRay.d) * Le) / (pdf * lightPdf); //initial

            if (!alpha.IsBlack()) {

                // Follow photon path through scene and record intersections
                PBRT_PHOTON_MAP_STARTED_RAY_PATH(&photonRay, &alpha);

                Intersection photonIsect;
                float vt0, vt1;
                //Photon Shooting & Depositing
                if (volume->IntersectP(photonRay, &vt0, &vt1));
            
                Point start = photonRay(vt0);
                Point end = photonRay(vt1);
                localPhotonBeams.push_back(PhotonBeam(start, end, alpha, end-start, 0.05f));
                
                PBRT_PHOTON_MAP_FINISHED_RAY_PATH(&photonRay, &alpha);
            }
        
            arena.FreeAll();
        
            // Merge local photon data with data in _PhotonIntegrator_
            MutexLock lock(mutex);
            
            // Give up if we're not storing enough photons
            if (abortTasks)
                return;
            if (nshot > 500000 && (unsuccessful(integrator->nPhotonBeamsWanted,
                                                PhotonBeams.size(), blockSize))){
                Error("Unable to store enough volume photons.  Giving up.\n");
                PhotonBeams.erase(PhotonBeams.begin(), PhotonBeams.end());
                abortTasks = true;
                return;
            }
            
            //update progress
            progress.Update(localPhotonBeams.size());
            nshot += blockSize;
            
            //Merge Volume Photons into main
            if (!volumeDone) {
                integrator->nVolumePaths += blockSize;
                for (uint32_t i = 0; i < localPhotonBeams.size(); ++i)
                    PhotonBeams.push_back(localPhotonBeams[i]);
                localPhotonBeams.erase(localPhotonBeams.begin(), localPhotonBeams.end());
                if (PhotonBeams.size() >= integrator->nPhotonBeamsWanted)
                    volumeDone = true;
            }
            
        }
        
        // Exit task if enough photons have been found
        if (volumeDone)
            break;
        }
    
}

void PhotonBeamIntegrator::Preprocess(const Scene *scene,
                                  const Camera *camera, const Renderer *renderer) {
    if (scene->lights.size() == 0) return;
    // Declare shared variables for photon shooting
    Mutex *mutex = Mutex::Create();
    vector<PhotonBeam> PhotonBeams;
    bool abortTasks = false;
    
    PhotonBeams.reserve(nPhotonBeamsWanted);
    uint32_t nshot = 0;
    
    // Compute light power CDF for photon shooting
    Distribution1D *lightDistribution = ComputeLightSamplingCDF(scene);
    
    // Run parallel tasks for photon shooting
    ProgressReporter progress(nPhotonBeamsWanted, "Shooting photons");
    vector<Task *> PhotonBeamShootingTasks;
    int nTasks = NumSystemCores();
    for (int i = 0; i < nTasks; ++i)
        PhotonBeamShootingTasks.push_back(new PhotonBeamShootingTask(i, camera ? camera->shutterOpen : 0.f, *mutex, this, progress, abortTasks,
                                                                    nshot, PhotonBeams, lightDistribution, scene, renderer));
    
    EnqueueTasks(PhotonBeamShootingTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < PhotonBeamShootingTasks.size(); ++i)
        delete PhotonBeamShootingTasks[i];
    
    Mutex::Destroy(mutex);
    progress.Done();
}

PhotonBeamIntegrator::PhotonBeamIntegrator(int nbeams)
{
    nPhotonBeamsWanted = nbeams;
}


PhotonBeamIntegrator *CreatePhotonBeamIntegrator(const ParamSet &params) {
    int nbeams = params.FindOneInt("photonbeams", 20000);
    return new PhotonBeamIntegrator(nbeams);
}
