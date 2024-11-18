#pragma once

#include "rsurface_types.h"
#include "surface_flow.h"
#include "remeshing/dynamic_remesher.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "scene_file.h"
#include "geometrycentral/surface/meshio.h"

#include "energy/squared_error.h"

#include <mkl.h>
#include "optimized_bct.h"
#include "energy/helfrich_energy.h"
#include "energy/helfrich_energy_ADE.h"
#include "energy/tpe_multipole_0.h"
#include "energy/tpe_barnes_hut_0.h"
#include "implicit/simple_surfaces.h"
#include "marchingcubes/CIsoSurface.h"

#define EIGEN_NO_DEBUG

namespace biorsurfaces
{
    struct PriorityVertex
    {
        GCVertex vertex;
        double priority;
        Vector3 position;
    };

    class MainApp
    {
    public:
        static int specifiedNumThreads;
        static int defaultNumThreads;
        static MainApp *instance;
        MainApp(MeshPtr mesh_, GeomPtr geom_, SurfaceFlow *flow_, polyscope::SurfaceMesh *psMesh_, std::string meshName_);

        void CreateAndDestroyBVH();
        void TestHelfrich();
        void TestMultiply();
        void TestUpdate();
        void TestObstacle0();
        void TestBarnesHut0();
        void PlotGradients();
        void Scale2x();
        void TestNormalDeriv();
        void MeshImplicitSurface(ImplicitSurface *surface);

        void GetFalloffWindow(GCVertex v, double radius, std::vector<PriorityVertex> &verts);
        void HandlePicking();

        void TakeOptimizationStep(bool remeshAfter, bool showAreaRatios);
        void AddObstacle(std::string filename, double weight, bool recenter, bool asPointCloud);
        void AddPotential(scene::PotentialType pType, double weight, double targetValue);
        void AddImplicitBarrier(scene::ImplicitBarrierData &implicitBarrier);

        MeshPtr mesh;
        GeomPtr geom;
        GeomPtr geomOrig;
        UVDataPtr uvs;
        SurfaceFlow *flow;
        TPEKernel *kernel;
        TPEnergyAllPairs *referenceEnergy;
        
        polyscope::SurfaceMesh *psMesh;
        std::vector<polyscope::SurfaceMesh *> obstacles;
        std::string meshName;
        int stepLimit=5000;
        int realTimeLimit;
        GradientMethod methodChoice;

        inline void reregisterMesh()
        {
            psMesh = polyscope::registerSurfaceMesh(meshName, geom->inputVertexPositions, mesh->getFaceVertexList(), polyscopePermutations(*mesh));
        }

        void updateMeshPositions();

        BVHNode6D *vertBVH;
        bool normalizeView;
        double bh_theta;
        remeshing::DynamicRemesher remesher;
        int numSteps;
        bool logPerformance;
        bool logEnergy;
        bool topoChange;
        int topoChangeSteps;
        long timeSpentSoFar;
        scene::SceneData sceneData;
        bool exitWhenDone;
        double totalObstacleVolume;
        bool ADE = false;
        std::vector<float> helfrichEnergyHistory; // Store Helfrich energy values
        std::vector<float> iterationHistory;         // Store iteration numbers
    private:
        int implicitCount = 0;
        GCVertex pickedVertex;
        std::vector<PriorityVertex> dragVertices;

        double pickDepth;
        void logPerformanceLine();
        void logEnergyFile();
        bool pickNearbyVertex(GCVertex &out);
        SquaredError *vertexPotential;
        bool ctrlMouseDown;
        Vector3 initialPickedPosition;
        bool hasPickedVertex;

    };
} // namespace biorsurfaces
