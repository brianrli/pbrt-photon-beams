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
    Point ipoint;
    Vector dir;
    float angle;
    float tcb; //photon beam distance
    float tbc; //ray distance
    float zmin;
    float zmax;
    float ri;
    int beamid;
    
    bool dup;
};

struct PhotonBeam {
    
    //cylinder
    PhotonBeam(const Point &st,
               const Point &en,
               const Spectrum &wt,//alpha
               const Vector &w, //direction
               float radius,
               int id)
    : start(st), end(en),alpha(wt), ri(radius),zmin(0),beamid(id),cone_radius(0),cone_height(0)
    {
        dir = Normalize(w);
        
        //necessary transformations
        WorldToObject = LookAt(st, end, Cross(end-start,Vector(1,0,0)));
        ObjectToWorld = Inverse(WorldToObject);
        zmax = (start-end).Length();
        
        scattered = true;
    }
    
    //cone
    PhotonBeam(const Point &st,
               const Point &en,
               const Spectrum &wt,//alpha
               const Vector &w, //direction
               float sa,
               int id,
               bool flag)
    : start(en), end(st),alpha(wt), zmin(0),beamid(id) //trade start and end
    {
        dir = Normalize(-w);
        
        //necessary transformations
        WorldToObject = LookAt(start, end, Cross(end-start,Vector(1,0,0)));
        ObjectToWorld = Inverse(WorldToObject);
        zmax = (start-end).Length();
                
        cone_height = zmax-zmin;
        cone_radius = cone_height * tanf(sa/2.f);
        
        scattered = false;
    }
    
    //split constructor
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
        beamid = beam.beamid;
        
        cone_radius = beam.cone_radius;
        cone_height = beam.cone_height;
        
        //new transformations
        WorldToObject = beam.WorldToObject;
        ObjectToWorld = beam.ObjectToWorld;
        
        //scattered beam?
        scattered = beam.scattered;
    }
    
    bool Intersect2(const Ray &r, beamisect &isect) const {
        float phi;
        Point phit;
        // Transform _Ray_ to object space
        Ray ray;
        WorldToObject(r, &ray);
        
        float radius = cone_radius;
        float height = cone_height;
        
        
        float A, B, C, t0, t1;
        if(!scattered){
            // Compute quadratic cone coefficients
            float k = radius / height;
            k = k*k;
            A = ray.d.x * ray.d.x + ray.d.y * ray.d.y -
                k * ray.d.z * ray.d.z;
            B = 2 * (ray.d.x * ray.o.x + ray.d.y * ray.o.y -
                           k * ray.d.z * (ray.o.z-height) );
            C = ray.o.x * ray.o.x + ray.o.y * ray.o.y -
            k * (ray.o.z -height) * (ray.o.z-height);
        }
        else{
            // Compute quadratic cylinder coefficients
            A = ray.d.x*ray.d.x + ray.d.y*ray.d.y;
            B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y);
            C = ray.o.x*ray.o.x + ray.o.y*ray.o.y - ri*ri;
        }
        
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
        
        // Compute cone inverse mapping
        phit = ray(thit);
        phi = atan2f(phit.y, phit.x);
        if (phi < 0.) phi += 2.f*M_PI;
        
        // Test cone intersection against clipping parameters
        if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
            if (thit == t1) return false;
            thit = t1;
            if (t1 > ray.maxt) return false;
            // Compute cone inverse mapping
            phit = ray(thit);
            phi = atan2f(phit.y, phit.x);
            if (phi < 0.) phi += 2.f*M_PI;
            if (phit.z < zmin || phit.z > zmax || phi > phiMax)
                return false;
        }
        
        //distance down ray to intersect point
        float raydist = thit;
        
        //local height+total height
        float beamdist = phit.z;
        //
        isect.ipoint = ObjectToWorld(phit);
        
        //compute angle
        Vector a = WorldToObject(r.d);
        Vector b = Vector(0,0,1);
        float angle = acosf(Dot(a,b)/(a.Length()+b.Length()));
        
        //angle, distance, power
        isect.angle = angle;
        isect.power = alpha;
        
        if(!scattered)
            isect.tcb = cone_height - beamdist;
        else
            isect.tcb = beamdist;
        
        if(!scattered){
            isect.ri = sqrtf(powf(phit.y,2) + powf(phit.x,2));
        }
        else
            isect.ri = ri;
        
        isect.tbc = raydist;
        isect.dir = dir;
        isect.beamid = beamid;
        isect.zmin = zmin;
        isect.zmax = zmax;
        isect.dup = false;
        return true;
    }

    
    BBox ObjectBound() const {
        Point p1;
        Point p2;
        if(!scattered){
            p1 = Point(-cone_radius, -cone_radius, zmin);
            p2 = Point( cone_radius,  cone_radius, zmax);
        }
        else{
            p1 = Point(-ri, -ri, zmin);
            p2 = Point( ri,  ri, zmax);
        }

        return ObjectToWorld(BBox(p1, p2));
    }
    
    /** Members **/
    Point start,end;
    Spectrum alpha;
    Vector dir;
    Transform ObjectToWorld, WorldToObject;
    float ri;
    int beamid;
    
    /** Cylinder **/
    float phiMax =2.f*M_PI;
    float zmin, zmax;
    
    /** Cone **/
    bool scattered;
    float cone_radius;
    float cone_height;
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
                    if (beams[node->primitivesOffset+i].Intersect2(ray,isect))
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
//    if(abs(sb.zmax - sb.zmin) < 1.f){
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

//use Henyey-Greenstein function
Vector VSampleHG(const Vector &w, float g,
                 float u1, float u2) {
	float costheta;
	if (fabsf(g) < 1e-3)
		costheta = 1.f - 2.f * u1;
	else
		costheta = -1.f / (2.f * g) *
        (1.f + g*g - ((1.f-g*g) * (1.f-g+2.f*g*u1)));
	float sintheta = sqrtf(max(0.f, 1.f-costheta*costheta));
	float phi = 2.f * M_PI * u2;
	Vector v1, v2;
	CoordinateSystem(w, &v1, &v2);
	return SphericalDirection(sintheta, costheta,
                              phi, v1, v2, w);
}

float VPhaseHG(const Vector &w, const Vector &wp, float g) {
    float costheta = Dot(w, wp);
    return 1.f / (4.f * M_PI) *
    (1.f - g*g) / powf(1.f + g*g - 2.f * g * costheta, 1.5f);
}

float SampleScattering(const Vector &wi, float u1, float u2, Vector &wo){
    wo = VSampleHG(-wi, .9995f, u1, u2);
    return VPhaseHG(-wi, wo, .9f);
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

    //unique id for each beam
    int beamid = integrator->nPhotonBeamsWanted * taskNum;

    while (true) {
//        const uint32_t blockSize = 500;
        const uint32_t blockSize = 25;
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

            float sa = 1/(pdf*integrator->nPhotonBeamsWanted);
            float aa = 2 * acosf(1 - sa/(2*M_PI));
            
            int paths=0;

            if (!alpha.IsBlack()) {

                // Follow photon path through scene and record intersections
                PBRT_PHOTON_MAP_STARTED_RAY_PATH(&photonRay, &alpha);

                Intersection photonIsect;
                float vt0, vt1;
                
                //Photon Shooting & Depositing
                int maxdepth = 3;
                int depth = 0;
                bool inVolume = true;
                float newradius;
                
                //intersect with group
                while (inVolume){
                    if (volume->IntersectP(photonRay, &vt0, &vt1))
                    {
                        if (scene->Intersect(photonRay, &photonIsect)) {
                            if(photonRay.maxt < vt1)
                                vt1 = photonRay.maxt;
                        }
                        
                        //deposit photonbeam
                        Point start = photonRay(vt0);
                        Point end = photonRay(vt1);
                        if(depth == 0){
                            //construct cone
                            split_beam(localPhotonBeams,PhotonBeam(start, end, alpha, Normalize(end-start), aa, beamid,true));
                        }
                        else
                            //construct cylinder
                            split_beam(localPhotonBeams,PhotonBeam(start, end, alpha, Normalize(end-start), aa, beamid));
                        paths++;
                        
                        //multiple scattering
                        float ts = (-logf(1-rng.RandomFloat())/volume->sigma_t(photonRay(vt0), photonRay.d, photonRay.time).y())+vt0;
                        
                        if(ts < vt1){

                            //scale down power
                            Spectrum st = volume->sigma_t(photonRay(ts), photonRay.d, photonRay.time);
                            Spectrum ss = volume->sigma_s(photonRay(ts), photonRay.d, photonRay.time);
                            alpha = alpha * (ts-vt0) * ss * Exp(-st * ts);
                            alpha *= (st/ss);
                            
                            //sample new direction
                            Vector newdirection;
                            float pdf = SampleScattering(photonRay.d, rng.RandomFloat(), rng.RandomFloat(),newdirection);
                            photonRay = RayDifferential(photonRay(ts),newdirection,0.f);
                            
//                            float cone_radius = (vt1-vt0) * tanf(aa/2.f);
//                            float coneheight = vt1-vt0;
                            if(depth==0){
                                aa = tanf(aa/2.f) * (vt1-ts);
                            }
                            
                            depth++;
                        }
                        else
                            break;
                        
                        //maximum depth check
                        if(depth >= maxdepth)
                            break;
                    }
                    else
                        break;
                }
                
                //avoid intersections of same beam
                beamid++;
                
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
                integrator->nVolumePaths += paths;
                for (uint32_t i = 0; i < localPhotonBeams.size(); ++i)
                    PhotonBeams.push_back(localPhotonBeams[i]);
                localPhotonBeams.erase(localPhotonBeams.begin(), localPhotonBeams.end());
                //based on how many beams
                if (integrator->nVolumePaths >= integrator->nPhotonBeamsWanted)
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
        PhotonBeamShootingTasks.push_back(new PhotonBeamShootingTask(i, camera ? camera->shutterOpen : 0.f, *mutex, this, progress, abortTasks,nshot, PhotonBeams, lightDistribution, scene, renderer,radius));
    
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

//biweight smoothing
float bi_kernel(float x){
    x = Clamp(x, 0.0f, 1.0f);
    float y = (15.f/16.f)*powf(1-powf(x,2),2.f);
    return y;
}

//Color change due to attenuation along the beam
Spectrum fb(Spectrum sigt, float v){
    return Exp(-sigt * v);
}

//color change due to attenuation towards the eye
Spectrum fe(Spectrum sigt, float u){
    return Exp(-sigt * u);
}

//shading depends on the viewing angle
Spectrum ff(Spectrum S, Spectrum p, float angle){
    return (S * (p/sinf(angle)));
}

//shading is influenced by the photon beam’s thickness
Spectrum ft(Spectrum power, float u){
    return (power * bi_kernel(u));
//    return power;
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
    intersections.reserve(10);
    BeamMap->Intersect(ray, intersections);
    
    float nints = intersections.size();
    vector<beamisect> copy = intersections;
    
    Spectrum L = Spectrum(0.f);
    Spectrum S;
    
    if(nints>0){
        for(int i = 0; i<intersections.size(); i++){
            for (int j = i+1; j< intersections.size(); j++){
                if (!intersections[i].dup
                    &&(intersections[i].beamid == intersections[j].beamid)){
                    if(intersections[j].tbc>intersections[i].tbc){
                        intersections[j].dup = true;
                    }
                    else{
                        intersections[i].dup = true;
                    }
                }
            }
        }
    }

    for (auto i = intersections.begin(); i != intersections.end(); i++) {
        
        //extinction coefficient
        if(!i->dup){
            Spectrum sigt = vr->sigma_t(i->ipoint,ray.d,ray.time);
            S = vr->sigma_s(i->ipoint,ray.d,ray.time);
            Spectrum P = vr->p(i->ipoint,i->dir,w,1e-4f);
            
    //        sigt = 0.8f; output 17
    //        sigt = 0.9f; //output 18
            L+=(ft(i->power,i->ri) * ff(S,P,i->angle) * fe(sigt,i->tbc) * fb(sigt,i->tcb));
        }
    }
    
    if(nints>0){
        L /= (ray.maxt-ray.mint);
    }
    
    return L;
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

PhotonBeamIntegrator::PhotonBeamIntegrator(int nbeams,float ssize,float beamradius)
{
    stepSize = ssize;
    nPhotonBeamsWanted = nbeams;
    radius = beamradius;
}


PhotonBeamIntegrator *CreatePhotonBeamIntegrator(const ParamSet &params) {
    int nbeams = params.FindOneInt("photonbeams", 20000);
    float stepSize  = params.FindOneFloat("stepsize", 4.f);
    float beamRadius  = params.FindOneFloat("beamradius", 0.001f);
    return new PhotonBeamIntegrator(nbeams,stepSize,beamRadius);
}
