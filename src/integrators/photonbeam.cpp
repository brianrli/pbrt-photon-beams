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
    PhotonBeam(const Point &st,
               const Point &en,
               const Spectrum &wt,//alpha
               const Vector &w, //direction
               float radius)
    : start(st), end(en),alpha(wt), ri(radius),zmin(0)
    {
        dir = Normalize(w);
        p = (st+en)/2.f;
        
        //necessary transformations
        WorldToObject = LookAt(st, end, Cross(end-start,Vector(1,0,0)));
        ObjectToWorld = Inverse(WorldToObject);
        zmax = (start-end).Length();
    }
    
    PhotonBeam(const Point &st,//need
               float zmi, //need
               PhotonBeam &beam)
    : start(st),zmin(zmi)
    {
        WorldToObject = beam.WorldToObject;
        ObjectToWorld = beam.ObjectToWorld;
        ri = beam.ri;
        dir = beam.dir;
        end = beam.end;
        zmax = beam.zmax;
        alpha = beam.alpha;
        p = (start+end)/2.f;
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
    
    BBox ObjectBound() const {
        Point p1 = Point(-ri, -ri, zmin);
        Point p2 = Point( ri,  ri, zmax);
        return ObjectToWorld(BBox(p1, p2));
    }
    

    /** Members **/
    Point start,end,p;
    Spectrum alpha;
    Vector dir;
    Transform ObjectToWorld, WorldToObject;
    float ri;
    BBox bound;
    
    /** Cylinder **/
    float phiMax =2.f*M_PI;
    float zmin, zmax;
};

void split_beam(vector<PhotonBeam> &beam,PhotonBeam sb)
{
    //subbeam has desired height
    if(abs(sb.zmax - sb.zmin) < 2.f){
        beam.push_back(sb);
    }
    //split beam in half
    else{
        float split_point = (sb.zmax - sb.zmin)/2.f;
        Point split = sb.start + (Normalize(sb.dir) * split_point);
        
        PhotonBeam newbeam = PhotonBeam(split,sb.zmin+split_point,sb);
        sb.end = split;
        sb.zmax = sb.zmin+split_point; //wrong
        
        split_beam(beam,sb);
        split_beam(beam,newbeam);
    }
}


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
                if (volume->IntersectP(photonRay, &vt0, &vt1))
                {
                    Point start = photonRay(vt0);
                    Point end = photonRay(vt1);
                    split_beam(localPhotonBeams,PhotonBeam(start, end, alpha, end-start, 0.05f));
                }
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
    
    //construct KdTree
    if (PhotonBeams.size() > 0)
        KdTree<PhotonBeam> *BeamMap = new KdTree<PhotonBeam>(PhotonBeams);
    
    
    
    bool flag = true;
}

PhotonBeamIntegrator::PhotonBeamIntegrator(int nbeams)
{
    nPhotonBeamsWanted = nbeams;
}


PhotonBeamIntegrator *CreatePhotonBeamIntegrator(const ParamSet &params) {
    int nbeams = params.FindOneInt("photonbeams", 20000);
    return new PhotonBeamIntegrator(nbeams);
}
