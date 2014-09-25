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

struct beamisect{
    beamisect()
    : power(Spectrum(0.f)),angle(0.f),tcb(0.f),tbc(0.f)
    {
    }
    
    Spectrum power;
    float angle;
    float tcb; //photon beam distance
    float tbc; //ray distance
};

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
        ri = beam.ri;
        dir = beam.dir;
        end = beam.end;
        zmax = beam.zmax;
        alpha = beam.alpha;

        //midpoint
        p = (start+end)/2.f;
        
        //new transformations
        WorldToObject = beam.WorldToObject;
        ObjectToWorld = beam.ObjectToWorld;
//        WorldToObject = LookAt(start, end, Cross(end-start,Vector(1,0,0)));
//        ObjectToWorld = Inverse(WorldToObject);
    }
    
    //computes intersect (bastardized cylinder code)
    bool Intersect(const Ray &r, beamisect &isect) const{
        float phi;
        Point phit;
        // Transform _Ray_ to object space
        Ray ray;
        WorldToObject(r, &ray);
        
        // Compute quadratic cylinder coefficients
        float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y;
        float B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y);
        float C = ray.o.x*ray.o.x + ray.o.y*ray.o.y - ri*ri;
        float t0, t1;
        
        // Solve quadratic equation for _t_ values
        if (!Quadratic(A, B, C, &t0, &t1))
            return false;
        
        // Compute intersection distance along ray
        if (t0 > ray.maxt || t1 < ray.mint)
            return false;
        float thit = t0;
        if (t0 < ray.mint) {
            thit = t1;
            if (thit > ray.maxt)
                return false;
        } //thit+t0
        
        float raydist = thit;
        
        // Compute cylinder hit point and $\phi$
        phit = ray(thit);
        
        //local height+total height
        float beamdist = phit.z + zmin;
        
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
        
        //compute angle
        Vector a = WorldToObject(dir);
        Vector b = Vector(0,0,1);
        float angle = acosf(Dot(a,b)/(a.Length()+b.Length()));
        
        //angle, distance, power
        isect.angle = angle;
        isect.power = alpha;
        isect.tcb = beamdist;
        isect.tbc = raydist;
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
    //BBox bound;
    
    /** Cylinder **/
    float phiMax =2.f*M_PI;
    float zmin, zmax;
};

// BBHAccel Local Declarations
struct BBHPrimitiveInfo {
    BBHPrimitiveInfo() { }
    BBHPrimitiveInfo(int pn, const BBox &b)
    : beamNumber(pn), bounds(b) {
        centroid = .5f * b.pMin + .5f * b.pMax;
    }
    int beamNumber;
    Point centroid;
    BBox bounds;
};


struct BBHBuildNode {
    // BBHBuildNode Public Methods
    BBHBuildNode() { children[0] = children[1] = NULL; }
    void InitLeaf(uint32_t first, uint32_t n, const BBox &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
    }
    void InitInterior(uint32_t axis, BBHBuildNode *c0, BBHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = Union(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
    }
    BBox bounds;
    BBHBuildNode *children[2];
    uint32_t splitAxis, firstPrimOffset, nPrimitives;
};

struct ComparePoints {
    ComparePoints(int d) { dim = d; }
    int dim;
    bool operator()(const BBHPrimitiveInfo &a,
                    const BBHPrimitiveInfo &b) const {
        return a.centroid[dim] < b.centroid[dim];
    }
};

struct CompareToBucket {
    CompareToBucket(int split, int num, int d, const BBox &b)
    : centroidBounds(b)
    { splitBucket = split; nBuckets = num; dim = d; }
    bool operator()(const BBHPrimitiveInfo &p) const;
    
    int splitBucket, nBuckets, dim;
    const BBox &centroidBounds;
};


bool CompareToBucket::operator()(const BBHPrimitiveInfo &p) const {
    int b = nBuckets * ((p.centroid[dim] - centroidBounds.pMin[dim]) /
                        (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
    if (b == nBuckets) b = nBuckets-1;
    Assert(b >= 0 && b < nBuckets);
    return b <= splitBucket;
}


struct LinearBBHNode {
    BBox bounds;
    union {
        uint32_t primitivesOffset;    // leaf
        uint32_t secondChildOffset;   // interior
    };
    
    uint8_t nPrimitives;  // 0 -> interior node
    uint8_t axis;         // interior node: xyz
    uint8_t pad[2];       // ensure 32 byte total size
};


static inline bool IntersectP(const BBox &bounds, const Ray &ray,
                              const Vector &invDir, const uint32_t dirIsNeg[3]) {
    // Check for ray intersection against $x$ and $y$ slabs
    float tmin =  (bounds[  dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tmax =  (bounds[1-dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tymin = (bounds[  dirIsNeg[1]].y - ray.o.y) * invDir.y;
    float tymax = (bounds[1-dirIsNeg[1]].y - ray.o.y) * invDir.y;
    if ((tmin > tymax) || (tymin > tmax))
        return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;
    
    // Check for ray intersection against $z$ slab
    float tzmin = (bounds[  dirIsNeg[2]].z - ray.o.z) * invDir.z;
    float tzmax = (bounds[1-dirIsNeg[2]].z - ray.o.z) * invDir.z;
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    return (tmin < ray.maxt) && (tmax > ray.mint);
}



// BBHAccel Method Definitions
BBHAccel::BBHAccel(const vector<PhotonBeam>&p, uint32_t mp)
{
    beams = p;
    maxBeamsInNode = min(255u, mp);
    
    if (beams.size() == 0) {
        nodes = NULL;
        return;
    }
    // Build BBH from _primitives_
    
    // Initialize _buildData_ array for primitives
    vector<BBHPrimitiveInfo> buildData;
    buildData.reserve(beams.size());
    for (uint32_t i = 0; i < beams.size(); ++i) {
        BBox bbox = beams[i].ObjectBound();
        buildData.push_back(BBHPrimitiveInfo(i, bbox));
    }
    
    // Recursively build BBH tree for primitives
    MemoryArena buildArena;
    uint32_t totalNodes = 0;
    vector<PhotonBeam> orderedBeams;
    orderedBeams.reserve(beams.size());
    BBHBuildNode *root = BBHAccel::recursiveBuild(buildArena, buildData, 0,
                                        beams.size(), &totalNodes,orderedBeams);
    beams.swap(orderedBeams);
    Info("BBH created with %d nodes for %d sub beams (%.2f MB)", totalNodes,
         (int)beams.size(), float(totalNodes * sizeof(LinearBBHNode))/(1024.f*1024.f));
    
    // Compute representation of depth-first traversal of BBH tree
    nodes = AllocAligned<LinearBBHNode>(totalNodes);
    for (uint32_t i = 0; i < totalNodes; ++i)
        new (&nodes[i]) LinearBBHNode;
    uint32_t offset = 0;
    flattenBBHTree(root, &offset);
    Assert(offset == totalNodes);
}


BBox BBHAccel::WorldBound() const {
    return nodes ? nodes[0].bounds : BBox();
}

BBHBuildNode *BBHAccel::recursiveBuild(MemoryArena &buildArena,
                                       vector<BBHPrimitiveInfo> &buildData, uint32_t start,
                                       uint32_t end, uint32_t *totalNodes,
                                       vector<PhotonBeam> &orderedBeams) {
    Assert(start != end);
    (*totalNodes)++; //might need this
    BBHBuildNode *node = buildArena.Alloc<BBHBuildNode>();
    
    // Compute bounds of all primitives in BBH node
    BBox bbox;
    for (uint32_t i = start; i < end; ++i)
        bbox = Union(bbox, buildData[i].bounds);
    uint32_t nPrimitives = end - start;
    if (nPrimitives == 1) {
        // Create leaf _BBHBuildNode_
        uint32_t firstPrimOffset = orderedBeams.size();
        for (uint32_t i = start; i < end; ++i) {
            uint32_t primNum = buildData[i].beamNumber;
            orderedBeams.push_back(beams[primNum]);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
    }
    else {
        // Compute bound of primitive centroids, choose split dimension _dim_
        BBox centroidBounds;
        for (uint32_t i = start; i < end; ++i)
            centroidBounds = Union(centroidBounds, buildData[i].centroid);
        int dim = centroidBounds.MaximumExtent();
        
        // Partition primitives into two sets and build children
        uint32_t mid = (start + end) / 2;
        if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
            node->InitInterior(dim,
                               recursiveBuild(buildArena, buildData, start, mid,
                                              totalNodes, orderedBeams),
                               recursiveBuild(buildArena, buildData, mid, end,
                                              totalNodes, orderedBeams));
            return node;
        }
        
        // Partition primitives based on _splitMethod_
        if (nPrimitives <= 4) {
            // Partition primitives into equally-sized subsets
            mid = (start + end) / 2;
            std::nth_element(&buildData[start], &buildData[mid],
                             &buildData[end-1]+1, ComparePoints(dim));
        }
        else {
            // Allocate _BucketInfo_ for SAH partition buckets
            const int nBuckets = 12;
            struct BucketInfo {
                BucketInfo() { count = 0; }
                int count;
                BBox bounds;
            };
            BucketInfo buckets[nBuckets];
            
            // Initialize _BucketInfo_ for SAH partition buckets
            for (uint32_t i = start; i < end; ++i) {
                int b = nBuckets *
                ((buildData[i].centroid[dim] - centroidBounds.pMin[dim]) /
                 (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
                if (b == nBuckets) b = nBuckets-1;
                Assert(b >= 0 && b < nBuckets);
                buckets[b].count++;
                buckets[b].bounds = Union(buckets[b].bounds, buildData[i].bounds);
            }
            
            // Compute costs for splitting after each bucket
            float cost[nBuckets-1];
            for (int i = 0; i < nBuckets-1; ++i) {
                BBox b0, b1;
                int count0 = 0, count1 = 0;
                for (int j = 0; j <= i; ++j) {
                    b0 = Union(b0, buckets[j].bounds);
                    count0 += buckets[j].count;
                }
                for (int j = i+1; j < nBuckets; ++j) {
                    b1 = Union(b1, buckets[j].bounds);
                    count1 += buckets[j].count;
                }
                cost[i] = .125f + (count0*b0.SurfaceArea() + count1*b1.SurfaceArea()) /
                bbox.SurfaceArea();
            }
            
            // Find bucket to split at that minimizes SAH metric
            float minCost = cost[0];
            uint32_t minCostSplit = 0;
            for (int i = 1; i < nBuckets-1; ++i) {
                if (cost[i] < minCost) {
                    minCost = cost[i];
                    minCostSplit = i;
                }
            }
            
            // Either create leaf or split primitives at selected SAH bucket
            if (nPrimitives > maxBeamsInNode ||
                minCost < nPrimitives) {
                BBHPrimitiveInfo *pmid = std::partition(&buildData[start],
                                                        &buildData[end-1]+1,
                                                        CompareToBucket(minCostSplit, nBuckets, dim, centroidBounds));
                mid = pmid - &buildData[0];
            }
            
            else {
                // Create leaf _BBHBuildNode_
                uint32_t firstPrimOffset = orderedBeams.size();
                for (uint32_t i = start; i < end; ++i) {
                    uint32_t primNum = buildData[i].beamNumber;
                    orderedBeams.push_back(beams[primNum]);
                }
                node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
                return node;
            }
        }
        
        node->InitInterior(dim,
                           recursiveBuild(buildArena, buildData, start, mid,
                                          totalNodes, orderedBeams),
                           recursiveBuild(buildArena, buildData, mid, end,
                                          totalNodes, orderedBeams));
    }
    return node;
}


uint32_t BBHAccel::flattenBBHTree(BBHBuildNode *node, uint32_t *offset) {
    LinearBBHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
    linearNode->nPrimitives = node->nPrimitives;
    uint32_t myOffset = (*offset)++;
    
    if (node->nPrimitives > 0) {
        Assert(!node->children[0] && !node->children[1]);
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    }
    else {
        // Creater interior flattened BBH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenBBHTree(node->children[0], offset);
        linearNode->secondChildOffset = flattenBBHTree(node->children[1],
                                                       offset);
    }
    return myOffset;
}

bool BBHAccel::Intersect(const Ray &ray, vector<beamisect> &intersections) const {
    if (!nodes) return false;
    bool hit = false;
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
    // Follow ray through BVH nodes to find primitive intersections
    uint32_t todoOffset = 0, nodeNum = 0;
    uint32_t todo[64];
    while (true) {
        const LinearBBHNode *node = &nodes[nodeNum];
        // Check ray against BVH node
        if (::IntersectP(node->bounds, ray, invDir, dirIsNeg)) {
            if (node->nPrimitives > 0) {
                // Intersect ray with primitives in leaf BVH node
                for (uint32_t i = 0; i < node->nPrimitives; ++i)
                {
                    beamisect isect;
                    if (beams[node->primitivesOffset+i].Intersect(ray,isect))
                    {
                        intersections.push_back(isect);
                        hit = true;
                    }
                }
                if (todoOffset == 0) break;
                nodeNum = todo[--todoOffset];
            }
            else {
                // Put far BVH node on _todo_ stack, advance to near node
                PBRT_BVH_INTERSECTION_TRAVERSED_INTERIOR_NODE(const_cast<LinearBVHNode *>(node));
                if (dirIsNeg[node->axis]) {
                    todo[todoOffset++] = nodeNum + 1;
                    nodeNum = node->secondChildOffset;
                }
                else {
                    todo[todoOffset++] = node->secondChildOffset;
                    nodeNum = nodeNum + 1;
                }
            }
        }
        else {
            if (todoOffset == 0) break;
            nodeNum = todo[--todoOffset];
        }
    }
    PBRT_BVH_INTERSECTION_FINISHED();
    return hit;
}


BBHAccel::~BBHAccel() {
    FreeAligned(nodes);
}

BBHAccel *CreateBBHAccelerator(const vector<PhotonBeam> &beams)
{
    return new BBHAccel(beams, 1);
}

void split_beam(vector<PhotonBeam> &beam,PhotonBeam sb)
{
    //subbeam has desired height
    if(abs(sb.zmax - sb.zmin) < 1.f){
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

                //intersect with group
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
        BeamMap = CreateBBHAccelerator(PhotonBeams);
}

Spectrum PhotonBeamIntegrator::Li(const Scene *scene,
            const Renderer *renderer,
            const RayDifferential &ray,
            const Sample *sample,
            RNG &rng,
            Spectrum *transmittance,
            MemoryArena &arena)const{
    if(!BeamMap){
		Error("Beam map is not initialized");
		exit(1);
	}
    VolumeRegion *vr = scene->volumeRegion;
    
    Vector w = -ray.d;
	vector<beamisect> intersections;
    intersections.reserve(100);
    BeamMap->Intersect(ray, intersections);
    
    float nints = intersections.size();
    
    float totalradius = 0;
    
    for (auto i = intersections.begin(); i != intersections.end(); i++) {
        return Spectrum(30.f);
        int stop = 1;
    }
    
    return Spectrum(0.f);
}

Spectrum PhotonBeamIntegrator::Transmittance(const Scene *scene,
                       const Renderer *,
                       const RayDifferential &ray,
                       const Sample *sample,
                       RNG &rng,
                       MemoryArena &arena)const
{
    if (!scene->volumeRegion) return Spectrum(1.f);
    float step, offset;
    if (sample) {
        step = stepSize;
        offset = sample->oneD[tauSampleOffset][0];
    }
    else {
        step = 4.f * stepSize;
        offset = rng.RandomFloat();
    }
    Spectrum tau = scene->volumeRegion->tau(ray, step, offset);
    return Exp(-tau);
}

void PhotonBeamIntegrator::RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene)
{
    tauSampleOffset = sample->Add1D(1);
	scatterSampleOffset = sample->Add1D(1);
}

PhotonBeamIntegrator::PhotonBeamIntegrator(int nbeams,float ssize)
{
    stepSize = ssize;
    nPhotonBeamsWanted = nbeams;
}


PhotonBeamIntegrator *CreatePhotonBeamIntegrator(const ParamSet &params) {
    int nbeams = params.FindOneInt("photonbeams", 20000);
    float stepSize  = params.FindOneFloat("stepsize", 4.f);
    return new PhotonBeamIntegrator(nbeams,stepSize);
}
