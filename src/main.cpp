//Note: Biorsurfaces is a software package for simulation of biological membranes. The majority of the code was derived from the c++ implementation of repulsive surfaces
//developed by Keenan Crane Lab at Carnegie Mellon. The package can be found at https://github.com/chrisyu-cs/repulsive-surfaces
//More information can be found at https://www.cs.cmu.edu/~kmcrane/Projects/RepulsiveSurfaces/index.html
//Main contributions of our project are the changes in customCallback() and addition of Helfrich energy and area difference elasticity terms and their analytical derivatives.
//Meghdad Razizadeh,St. Jude Children's Research Hospital, 2024


#include "main.h"
#include "main_picking.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "polyscope/point_cloud.h"
#include "../deps/polyscope/deps/args/args/args.hxx"
#include "imgui.h"
#include "surface_derivatives.h"
#include "energy/tpe_kernel.h"
#include "energy/all_energies.h"
#include "helpers.h"
#include <memory>
#include <Eigen/Sparse>
#include <omp.h>
#include "sobolev/all_constraints.h"
#include "sobolev/hs.h"
#include "sobolev/hs_iterative.h"
#include "sobolev/h1.h"
#include "spatial/convolution.h"
#include "spatial/convolution_kernel.h"
#include "surface_derivatives.h"
#include "obj_writer.h"
#include "dropdown_strings.h"
#include "energy/coulomb.h"
#include "energy/helfrich_energy.h"
#include "energy/helfrich_energy_ADE.h"
#include "bct_constructors.h"
#include "remeshing/remeshing.h"
#include "imgui.h"  // Include ImGui first
#include "line_search.h"
using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace biorsurfaces
{
    int MainApp::specifiedNumThreads;
    int MainApp::defaultNumThreads;
    MainApp *MainApp::instance = 0;

    MainApp::MainApp(MeshPtr mesh_, GeomPtr geom_, SurfaceFlow *flow_, polyscope::SurfaceMesh *psMesh_, std::string meshName_)
        : mesh(std::move(mesh_)), geom(std::move(geom_)), geomOrig(geom->copy()), remesher(mesh, geom, geomOrig)
    {
        flow = flow_;
        psMesh = psMesh_;
        meshName = meshName_;
        vertBVH = 0;
        vertexPotential = 0;
        ctrlMouseDown = false;
        hasPickedVertex = false;
        numSteps = 0;
        methodChoice = GradientMethod::HsProjectedIterative;
        timeSpentSoFar = 0;
        realTimeLimit = 0;
        logEnergy = true;
        topoChange = true;
        topoChangeSteps = 5;
        referenceEnergy = 0;
        exitWhenDone = false;
        totalObstacleVolume = 0;
    }

    void MainApp::logEnergyFile()
    {
        std::ofstream outputFile("log_Energy.txt", std::ios::app); // Open output file in append mode.
        outputFile << biorsurfaces::MainApp::instance->flow->energies[1]->Value() / biorsurfaces::MainApp::instance->flow->energies[1]->GetWeight() << " "
                   << biorsurfaces::MainApp::instance->flow->energies[0]->Value() / biorsurfaces::MainApp::instance->flow->energies[0]->GetWeight() << " "
                   << std::endl;
    }

    void MainApp::TakeOptimizationStep(bool remeshAfter, bool showAreaRatios)
    {
        ptic("MainApp::TakeOptimizationStep");

        if (logEnergy)
        {
            logEnergyFile();
        }

        long beforeStep = currentTimeMilliseconds();

        ptic("Switch");
        switch (methodChoice)
        {
        case GradientMethod::HsProjected:
            flow->StepProjectedGradient();
            break;
        case GradientMethod::HsProjectedIterative:
            flow->StepProjectedGradientIterative();
            break;
        case GradientMethod::HsExactProjected:
            flow->StepProjectedGradientExact();
            break;
        case GradientMethod::H1Projected:
            flow->StepH1ProjGrad();
            break;
        case GradientMethod::L2Unconstrained:
            flow->StepL2Unconstrained();
            break;
        case GradientMethod::L2Projected:
            flow->StepL2Projected();
            break;
        case GradientMethod::AQP:
        {
            double kappa = 100;
            flow->StepAQP(1 / kappa);
        }
        break;
        case GradientMethod::H1_LBFGS:
            flow->StepH1LBFGS();
            break;
        case GradientMethod::BQN_LBFGS:
            flow->StepBQN();
            break;
        case GradientMethod::H2Projected:
        case GradientMethod::Helfrich:
            flow->StepH2Projected();
            break;
        default:
            throw std::runtime_error("Unknown gradient method type.");
        }
        ptoc("Switch");

        if (remeshAfter)
        {
            bool doCollapse = (topoChange && (numSteps % topoChangeSteps == 0));
            std::cout << "Applying remeshing..." << std::endl;
            flow->verticesMutated = remesher.Remesh(5, doCollapse);
            if (flow->verticesMutated)
            {
                std::cout << "Vertices were mutated this step -- memory vectors are now invalid." << std::endl;
            }
            else
            {
                std::cout << "Vertices were not mutated this step." << std::endl;
            }
            ptic("mesh->compress()");
            mesh->compress();
            ptoc("mesh->compress()");
            ptic("MainApp::instance->reregisterMesh();");
            MainApp::instance->reregisterMesh();
            ptoc("MainApp::instance->reregisterMesh();");
        }
        else
        {
            flow->verticesMutated = false;
            MainApp::instance->updateMeshPositions();
        }
        long afterStep = currentTimeMilliseconds();
        long timeForStep = afterStep - beforeStep;
        timeSpentSoFar += timeForStep;
        numSteps++;
        std::cout << "  Mesh total volume = " << totalVolume(geom, mesh) << std::endl;
        std::cout << "  Mesh total area = " << totalArea(geom, mesh) << std::endl;

        if (showAreaRatios)
        {
            VertexData<double> areaRatio(*mesh);
            for (Vertex v : mesh->vertices())
            {
                areaRatio[v] = geomOrig->vertexDualArea(v) / geom->vertexDualArea(v);
            }

            psMesh->addVertexScalarQuantity("Area ratios", areaRatio);
        }
        float currentHelfrichEnergy = flow->energies[1]->Value();
        helfrichEnergyHistory.push_back(currentHelfrichEnergy);
        iterationHistory.push_back(numSteps); 

        ptoc("MainApp::TakeOptimizationStep");
    }

    void MainApp::updateMeshPositions()
    {
        if (normalizeView)
        {
            double scale = 0;
            for (GCVertex v : mesh->vertices())
            {
                scale = fmax(scale, norm(geom->inputVertexPositions[v]));
            }
            std::vector<Vector3> scaled(mesh->nVertices());
            VertexIndices inds = mesh->getVertexIndices();
            for (GCVertex v : mesh->vertices())
            {
                scaled[inds[v]] = geom->inputVertexPositions[v] / scale;
            }
            psMesh->updateVertexPositions(scaled);
        }
        else
        {
            psMesh->updateVertexPositions(geom->inputVertexPositions);
        }
        polyscope::requestRedraw();
    }

    void PlotMatrix(Eigen::MatrixXd &mat, polyscope::SurfaceMesh *psMesh, std::string name)
    {
        std::vector<Vector3> vecs;
        for (int i = 0; i < mat.rows(); i++)
        {
            Vector3 row_i = GetRow(mat, i);
            vecs.push_back(row_i);
        }
        psMesh->addVertexVectorQuantity(name, vecs);
    }

    void PlotVector(Eigen::VectorXd &vec, int nVerts, polyscope::SurfaceMesh *psMesh, std::string name)
    {
        Eigen::MatrixXd M;
        M.setZero(nVerts, 3);
        MatrixUtils::ColumnIntoMatrix(vec, M);
        PlotMatrix(M, psMesh, name);
    }

    void MainApp::PlotGradients()
    {
        Eigen::MatrixXd l2Diff, hsGrad, hsGradExact,helfGrad, TPEGrad;
        l2Diff.setZero(mesh->nVertices(), 3);
        hsGrad.setZero(mesh->nVertices(), 3);
        hsGradExact.setZero(mesh->nVertices(), 3);
        TPEGrad.setZero(mesh->nVertices(), 3);
        helfGrad.setZero(mesh->nVertices(), 3);
        flow->UpdateEnergies();

        std::cout << "Assembling L2 differential..." << std::endl;
        long diffTimeStart = currentTimeMilliseconds();
        flow->AssembleGradients(l2Diff);
        long diffTimeEnd = currentTimeMilliseconds();
        std::cout << "Differential took " << (diffTimeEnd - diffTimeStart) << " ms" << std::endl;

        // std::unique_ptr<Hs::HsMetric> hs = flow->GetHsMetric();

        // std::cout << "Inverting \"sparse\" metric..." << std::endl;
        // long sparseTimeStart = currentTimeMilliseconds();
        // hs->InvertMetricMat(l2Diff, hsGrad);
        // long sparseTimeEnd = currentTimeMilliseconds();
        // std::cout << "Sparse metric took " << (sparseTimeEnd - sparseTimeStart) << " ms" << std::endl;

        // std::cout << "Inverting dense metric..." << std::endl;
        // long timeStart = currentTimeMilliseconds();
        // std::vector<ConstraintPack> empty;
        // // hs->ProjectGradientExact(l2Diff, hsGradExact, empty);
        // hsGradExact = hsGrad;
        // long timeEnd = currentTimeMilliseconds();
        // std::cout << "Dense metric took " << (timeEnd - timeStart) << " ms" << std::endl;

        std::cout << "Assembling Membrane Energy Gradient" << std::endl;
        long diffTimeStart1 = currentTimeMilliseconds();
        biorsurfaces::MainApp::instance->flow->energies[1]->Differential(helfGrad);      
        long diffTimeEnd1 = currentTimeMilliseconds();
        std::cout << "Differential took " << (diffTimeEnd1 - diffTimeStart1) << " ms" << std::endl;

        std::cout << "Assembling TPE Gradient" << std::endl;
        long diffTimeStart2 = currentTimeMilliseconds();
        biorsurfaces::MainApp::instance->flow->energies[0]->Differential(TPEGrad);      
        long diffTimeEnd2 = currentTimeMilliseconds();
        std::cout << "Differential took " << (diffTimeEnd2 - diffTimeStart2) << " ms" << std::endl;



        PlotMatrix(l2Diff, psMesh, "L2 differential");
        // PlotMatrix(hsGrad, psMesh, "Hs sparse gradient");
        // PlotMatrix(hsGradExact, psMesh, "Hs dense gradient");
        PlotMatrix(TPEGrad, psMesh, "Tangent Point Energy gradient");
        PlotMatrix(helfGrad, psMesh, "Membrane Energy gradient");

    }

    bool MainApp::pickNearbyVertex(GCVertex &out)
    {
        using namespace polyscope;
        Vector2 screenPos = getMouseScreenPos();

        std::pair<Structure *, size_t> pickVal =
            pick::evaluatePickQuery(screenPos.x, screenPos.y);

        GCVertex pickedVert;
        GCFace pickedFace;
        GCEdge pickedEdge;
        GCHalfedge pickedHalfedge;

        glm::mat4 view = polyscope::view::getCameraViewMatrix();
        glm::mat4 proj = polyscope::view::getCameraPerspectiveMatrix();
        glm::mat4 viewProj = proj * view;

        polyscope::SurfaceMesh *asMesh = dynamic_cast<polyscope::SurfaceMesh *>(pickVal.first);

        if (tryGetPickedVertex(asMesh, pickVal.second, mesh, pickedVert))
        {
            out = pickedVert;
            return true;
        }
        else if (tryGetPickedFace(asMesh, pickVal.second, mesh, pickedFace))
        {
            out = nearestVertexToScreenPos(screenPos, geom, viewProj, pickedFace);
            return true;
        }
        else if (tryGetPickedEdge(asMesh, pickVal.second, mesh, pickedEdge))
        {
            out = nearestVertexToScreenPos(screenPos, geom, viewProj, pickedEdge);
            return true;
        }
        else if (tryGetPickedHalfedge(asMesh, pickVal.second, mesh, pickedHalfedge))
        {
            out = nearestVertexToScreenPos(screenPos, geom, viewProj, pickedHalfedge);
            return true;
        }
        else
        {
            std::cout << "No valid element was picked (index " << pickVal.second << ")" << std::endl;
            return false;
        }
    }

    class PVCompare
    {
    public:
        bool operator()(PriorityVertex v1, PriorityVertex v2)
        {
            return (v1.priority > v2.priority);
        }
    };

    double gaussian(double radius, double dist)
    {
        double radterm = dist / radius;
        double epow = exp(-0.5 * radterm * radterm);
        return epow;
    }

    void MainApp::GetFalloffWindow(GCVertex v, double radius, std::vector<PriorityVertex> &verts)
    {
        // Do a simple Dijkstra search on edges
        VertexData<bool> seen(*mesh, false);
        std::priority_queue<PriorityVertex, std::vector<PriorityVertex>, PVCompare> queue;
        queue.push(PriorityVertex{v, 0, geom->inputVertexPositions[v]});

        while (!queue.empty())
        {
            PriorityVertex next = queue.top();
            queue.pop();

            if (next.priority > radius)
            {
                break;
            }
            else if (seen[next.vertex])
            {
                continue;
            }
            else
            {
                // Mark the next vertex as seen
                seen[next.vertex] = true;
                // Compute the weight
                double weight = gaussian(radius / 3, next.priority);
                verts.push_back(PriorityVertex{next.vertex, weight, geom->inputVertexPositions[next.vertex]});

                // Enqueue all neighbors
                for (GCVertex neighbor : next.vertex.adjacentVertices())
                {
                    if (seen[neighbor])
                    {
                        continue;
                    }
                    // Add the next edge distance
                    Vector3 p1 = geom->inputVertexPositions[next.vertex];
                    Vector3 p2 = geom->inputVertexPositions[neighbor];
                    double neighborDist = next.priority + norm(p1 - p2);

                    queue.push(PriorityVertex{neighbor, neighborDist, geom->inputVertexPositions[neighbor]});
                }
            }
        }

        std::cout << "Got " << verts.size() << " vertices" << std::endl;
    }

    void MainApp::HandlePicking()
    {
        using namespace polyscope;

        auto io = ImGui::GetIO();
        glm::mat4 view = polyscope::view::getCameraViewMatrix();
        glm::mat4 proj = polyscope::view::getCameraPerspectiveMatrix();
        glm::mat4 viewProj = proj * view;

        if (io.KeyCtrl && io.MouseDown[0])
        {
            if (!ctrlMouseDown)
            {
                if (pickNearbyVertex(pickedVertex))
                {
                    hasPickedVertex = true;
                    GetFalloffWindow(pickedVertex, 0.5, dragVertices);

                    Vector3 screen = projectToScreenCoords3(geom->inputVertexPositions[pickedVertex], viewProj);
                    pickDepth = screen.z;

                    Vector3 unprojected = unprojectFromScreenCoords3(Vector2{screen.x, screen.y}, pickDepth, viewProj);
                    initialPickedPosition = geom->inputVertexPositions[pickedVertex];
                }
                ctrlMouseDown = true;
            }
            else
            {
                if (hasPickedVertex)
                {
                    Vector2 mousePos = getMouseScreenPos();
                    Vector3 unprojected = unprojectFromScreenCoords3(mousePos, pickDepth, viewProj);
                    Vector3 displacement = unprojected - initialPickedPosition;

                    for (PriorityVertex &v : dragVertices)
                    {
                        Vector3 newPos = v.position + v.priority * displacement;
                        geom->inputVertexPositions[v.vertex] = newPos;
                    }

                    flow->ResetAllConstraints();
                    // INACTIVATED BY MEGHDAD-Jan 3, 2023

                    flow->ResetAllPotentials();

                    if (vertexPotential)
                    {
                        for (PriorityVertex &v : dragVertices)
                        {
                            vertexPotential->ChangeVertexTarget(v.vertex, geom->inputVertexPositions[v.vertex]);
                        }
                    }
                    updateMeshPositions();
                }
            }
        }
        else
        {
            if (ctrlMouseDown)
            {
                ctrlMouseDown = false;
                hasPickedVertex = false;
                dragVertices.clear();
                geom->inputVertexPositions[pickedVertex] = initialPickedPosition;
                updateMeshPositions();
            }
        }
    }

    void MainApp::CreateAndDestroyBVH()
    {
        OptimizedClusterTree *bvh = CreateOptimizedBVH(mesh, geom);
        std::cout << "Created BVH" << std::endl;
        delete bvh;
        std::cout << "Deleted BVH" << std::endl;
    }

    void MainApp::Scale2x()
    {
        for (GCVertex v : mesh->vertices())
        {
            geom->inputVertexPositions[v] = 2 * geom->inputVertexPositions[v];
        }
    }

    class VectorInit
    {
    public:
        static void Init(Vector3 &data, BVHNode6D *node)
        {
            data = Vector3{1, 2, 3};
        }
    };

    void MainApp::AddObstacle(std::string filename, double weight, bool recenter, bool asPointCloud)
    {
        std::unique_ptr<surface::SurfaceMesh> obstacleMesh;
        GeomUPtr obstacleGeometry;
        // Load mesh
        std::tie(obstacleMesh, obstacleGeometry) = readNonManifoldMesh(filename);

        obstacleGeometry->requireVertexDualAreas();
        obstacleGeometry->requireVertexNormals();

        if (recenter)
        {
            Vector3 obstacleCenter = meshBarycenter(obstacleGeometry, obstacleMesh);
            std::cout << "Recentering obstacle " << filename << " (offset " << obstacleCenter << ")" << std::endl;
            for (GCVertex v : obstacleMesh->vertices())
            {
                obstacleGeometry->inputVertexPositions[v] = obstacleGeometry->inputVertexPositions[v] - obstacleCenter;
            }
        }

        std::string mesh_name = polyscope::guessNiceNameFromPath(filename);

        if (asPointCloud)
        {
            polyscope::PointCloud *pointCloud = polyscope::registerPointCloud(mesh_name, obstacleGeometry->inputVertexPositions);
        }

        else
        {
            polyscope::SurfaceMesh *psMesh = polyscope::registerSurfaceMesh(mesh_name, obstacleGeometry->inputVertexPositions,
                                                                            obstacleMesh->getFaceVertexList(), polyscopePermutations(*obstacleMesh));
        }

        std::unique_ptr<surface::SurfaceMesh> sharedObsMesh = std::move(obstacleMesh);
        GeomPtr sharedObsGeom = std::move(obstacleGeometry);

        SurfaceEnergy *obstacleEnergy = 0;

        if (asPointCloud)
        {
            size_t nVerts = sharedObsMesh->nVertices();

            Eigen::VectorXd wts;
            wts.setOnes(nVerts);

            Eigen::MatrixXd pos;
            pos.setZero(nVerts, 3);

            for (size_t i = 0; i < nVerts; i++)
            {
                Vector3 v = sharedObsGeom->inputVertexPositions[i];
                MatrixUtils::SetRowFromVector3(pos, i, v);
            }

            obstacleEnergy = new TPPointCloudObstacleBarnesHut0(mesh, geom, flow->BaseEnergy(), wts, pos,
                                                                kernel->alpha, kernel->beta, bh_theta, weight);
        }

        else
        {
            obstacleEnergy = new TPObstacleBarnesHut0(mesh, geom, flow->BaseEnergy(), sharedObsMesh, sharedObsGeom,
                                                      kernel->alpha, kernel->beta, bh_theta, weight);
        }

        flow->AddObstacleEnergy(obstacleEnergy);
        std::cout << "Added " << filename << " as obstacle with weight " << weight << std::endl;

        totalObstacleVolume += totalVolume(sharedObsGeom, sharedObsMesh);
    }

    void MainApp::AddImplicitBarrier(scene::ImplicitBarrierData &barrierData)
    {
        ImplicitSurface *implSurface;
        // Create the requested implicit surface
        switch (barrierData.type)
        {
        case scene::ImplicitType::Plane:
        {
            Vector3 point{barrierData.parameters[0], barrierData.parameters[1], barrierData.parameters[2]};
            Vector3 normal{barrierData.parameters[3], barrierData.parameters[4], barrierData.parameters[5]};
            std::cout << "Constructing implicit plane at point " << point << " with normal " << normal << std::endl;
            implSurface = new FlatPlane(point, normal);
        }
        break;
        case scene::ImplicitType::Torus:
        {
            double major = barrierData.parameters[0];
            double minor = barrierData.parameters[1];
            Vector3 center{barrierData.parameters[2], barrierData.parameters[3], barrierData.parameters[4]};
            std::cout << "Constructing implicit torus with major radius " << major << ", minor radius " << minor << ", center " << center << std::endl;
            implSurface = new ImplicitTorus(major, minor, center);
        }
        break;
        case scene::ImplicitType::Sphere:
        {
            double radius = barrierData.parameters[0];
            Vector3 center{barrierData.parameters[1], barrierData.parameters[2], barrierData.parameters[3]};
            std::cout << "Constructing implicit sphere with radius " << radius << ", center " << center << std::endl;
            implSurface = new ImplicitSphere(radius, center);
        }
        break;
        case scene::ImplicitType::Cylinder:
        {
            double radius = barrierData.parameters[0];
            Vector3 center{barrierData.parameters[1], barrierData.parameters[2], barrierData.parameters[3]};
            Vector3 axis{barrierData.parameters[4], barrierData.parameters[5], barrierData.parameters[6]};
            std::cout << "Constructing implicit cylinder with radius " << radius << ", center " << center << ", axis " << axis << std::endl;
            implSurface = new ImplicitCylinder(radius, center, axis);
        }
        break;
        default:
        {
            throw std::runtime_error("Unimplemented implicit surface type.");
        }
        break;
        }

        // Mesh the 0 isosurface so we can see the implicit surface
        MainApp::instance->MeshImplicitSurface(implSurface);

        // Use the implicit surface to setup the energy
        std::unique_ptr<ImplicitSurface> implUnique(implSurface);
        if (barrierData.repel)
        {
            std::cout << "Using implicit surface as obstacle, with power " << barrierData.power << " and weight " << barrierData.weight << std::endl;
            ImplicitObstacle *obstacle = new ImplicitObstacle(mesh, geom, std::move(implUnique), barrierData.power, barrierData.weight);
            flow->AddAdditionalEnergy(obstacle);
        }
        else
        {
            std::cout << "Using implicit surface as attractor, with power " << barrierData.power << " and weight " << barrierData.weight << std::endl;
            ImplicitAttractor *attractor = new ImplicitAttractor(mesh, geom, std::move(implUnique), uvs, barrierData.power, barrierData.weight);
            flow->AddAdditionalEnergy(attractor);
        }
    }

    void MainApp::AddPotential(scene::PotentialType pType, double weight, double targetValue)
    {
        switch (pType)
        {
        case scene::PotentialType::SquaredError:
        {
            SquaredError *errorPotential = new SquaredError(mesh, geom, weight);
            vertexPotential = errorPotential;
            flow->AddAdditionalEnergy(errorPotential);
            remesher.KeepVertexDataUpdated(&errorPotential->originalPositions);
            break;
        }
        case scene::PotentialType::Area:
        {
            TotalAreaPotential *areaPotential = new TotalAreaPotential(mesh, geom, weight);
            flow->AddAdditionalEnergy(areaPotential);
            break;
        }
        case scene::PotentialType::Volume:
        {
            TotalVolumePotential *volumePotential = new TotalVolumePotential(mesh, geom, weight);
            flow->AddAdditionalEnergy(volumePotential);
            break;
        }
        case scene::PotentialType::BoundaryLength:
        {
            BoundaryLengthPenalty *errorPotential = new BoundaryLengthPenalty(mesh, geom, weight, targetValue);
            flow->AddAdditionalEnergy(errorPotential);
            break;
        }
        case scene::PotentialType::BoundaryCurvature:
        {
            BoundaryCurvaturePenalty *errorPotential = new BoundaryCurvaturePenalty(mesh, geom, weight);
            flow->AddAdditionalEnergy(errorPotential);
            break;
        }
        case scene::PotentialType::SoftAreaConstraint:
        {
            SoftAreaConstraint *softArea = new SoftAreaConstraint(mesh, geom, weight);
            flow->AddAdditionalEnergy(softArea);
            break;
        }
        case scene::PotentialType::SoftVolumeConstraint:
        {
            SoftVolumeConstraint *softVol = new SoftVolumeConstraint(mesh, geom, weight);
            flow->AddAdditionalEnergy(softVol);
            break;
        }
        case scene::PotentialType::Helfrich:
        {
            HelfrichEnergy *helfrich = new HelfrichEnergy(mesh, geom, weight);
            flow->AddAdditionalEnergy(helfrich);
            break;
        }

        case scene::PotentialType::HelfrichADE:
        {
            HelfrichEnergyADE *helfrichADE = new HelfrichEnergyADE(mesh, geom, weight);
            flow->AddAdditionalEnergy(helfrichADE);
            ADE = true;
            break;
        }
        case scene::PotentialType::SoftReducedVolumeConstraint:
        {
            SoftReducedVolumeConstraint *softRV = new SoftReducedVolumeConstraint(mesh, geom, weight, targetValue);
            flow->AddAdditionalEnergy(softRV);
            break;
        }

        case scene::PotentialType::SoftTotalMeanCurvatureConstraint:
        {
            SoftTotalMeanCurvatureConstraint *softTMC = new SoftTotalMeanCurvatureConstraint(mesh, geom, weight, targetValue);
            flow->AddAdditionalEnergy(softTMC);
            break;
        }

        default:
        {
            std::cout << "Unknown potential type." << std::endl;
            break;
        }
        }
    }

    void MainApp::MeshImplicitSurface(ImplicitSurface *surface)
    {
        CIsoSurface<double> *iso = new CIsoSurface<double>();

        std::cout << "Meshing the supplied implicit surface using marching cubes..." << std::endl;

        const int numCells = 50;
        Vector3 center = surface->BoundingCenter();
        double diameter = surface->BoundingDiameter();
        double cellSize = diameter / numCells;
        double radius = diameter / 2;

        Vector3 lowerCorner = center - Vector3{radius, radius, radius};

        int numCorners = numCells + 1;

        double field[numCorners * numCorners * numCorners];

        int nSlice = numCorners * numCorners;
        int nRow = numCorners;

        for (int x = 0; x < numCorners; x++)
        {
            for (int y = 0; y < numCorners; y++)
            {
                for (int z = 0; z < numCorners; z++)
                {
                    Vector3 samplePt = lowerCorner + Vector3{(double)x, (double)y, (double)z} * cellSize;
                    double value = surface->SignedDistance(samplePt);
                    field[nSlice * z + nRow * y + x] = value;
                }
            }
        }

        iso->GenerateSurface(field, 0, numCells, numCells, numCells, cellSize, cellSize, cellSize);

        std::vector<glm::vec3> nodes;
        std::vector<std::array<size_t, 3>> triangles;

        int nVerts = iso->m_nVertices;

        for (int i = 0; i < nVerts; i++)
        {
            double x = iso->m_ppt3dVertices[i][0];
            double y = iso->m_ppt3dVertices[i][1];
            double z = iso->m_ppt3dVertices[i][2];

            Vector3 p = lowerCorner + Vector3{x, y, z};
            nodes.push_back(glm::vec3{p.x, p.y, p.z});
        }

        int nTris = iso->m_nTriangles;

        for (int i = 0; i < nTris; i++)
        {
            int i1 = iso->m_piTriangleIndices[3 * i];
            int i2 = iso->m_piTriangleIndices[3 * i + 1];
            int i3 = iso->m_piTriangleIndices[3 * i + 2];

            triangles.push_back({(size_t)i1, (size_t)i2, (size_t)i3});
        }

        implicitCount++;
        polyscope::registerSurfaceMesh("implicitSurface" + std::to_string(implicitCount), nodes, triangles);
        delete iso;
    }
} // namespace biorsurfaces

// UI parameters
bool run = false;
bool takeScreenshots = false;
bool saveOBJs = false;
uint screenshotNum = 0;
uint objNum = 0;
bool uiNormalizeView = false;
bool remesh = true;
bool changeTopo = true;
bool areaRatios = false;
bool ADE = true;
bool limitStep = false;
float maxStep = -1.;
float regMC = 1.2;
float regWMC = 1e7;
float regWW = 5000;
float changeStep = 0.025;
float regalpha = 2/ geometrycentral::PI;
float regCons = 1.05;
float regConsArea = 0.75;
int partIndex = 4475;

enum class ModelType {
    BC,
    ADE
};
ModelType currentModel = ModelType::ADE;

void saveScreenshot(uint i)
{
    char buffer[5];
    std::snprintf(buffer, sizeof(buffer), "%04d", i);
    std::string fname = "frames/frame" + std::string(buffer) + ".png";
    polyscope::screenshot(fname, false);
    std::cout << "Saved screenshot to " << fname << std::endl;
}

void saveOBJ(biorsurfaces::MeshPtr mesh, biorsurfaces::GeomPtr geom, biorsurfaces::GeomPtr geomOrig, uint i)
{

    char buffer[5];
    std::snprintf(buffer, sizeof(buffer), "%04d", i);
    std::string fname = "mesh/frame" + std::string(buffer) + ".obj";
    biorsurfaces::writeMeshToOBJ(mesh, geom, geomOrig, areaRatios, fname);
    std::cout << "Saved OBJ frame to " << fname << std::endl;
}

template <typename ItemType>
void selectFromDropdown(std::string label, const ItemType choices[], size_t nChoices, ItemType &store)
{
    using namespace biorsurfaces;

    // Dropdown menu for list of remeshing mode settings
    if (ImGui::BeginCombo(label.c_str(), StringOfMode(store).c_str()))
    {
        for (size_t i = 0; i < nChoices; i++)
        {
            bool is_selected = (store == choices[i]);
            if (ImGui::Selectable(StringOfMode(choices[i]).c_str(), is_selected))
                store = choices[i];
            if (is_selected)
                ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h


void customCallback()
{
    using namespace biorsurfaces;
    const int INDENT = 10;
    const int ITEM_WIDTH = 160;
    ImGui::Text("Flow control");
    ImGui::BeginGroup();
    ImGui::Indent(INDENT);
    ImGui::PushItemWidth(ITEM_WIDTH);
    ImGui::Checkbox("Run flow", &run);
    ImGui::SameLine(ITEM_WIDTH, 2 * INDENT);
    ImGui::Checkbox("Normalize view", &uiNormalizeView);
    ImGui::Checkbox("Take screenshots", &takeScreenshots);
    ImGui::SameLine(ITEM_WIDTH, 2 * INDENT);
    if ((takeScreenshots && screenshotNum == 0) || ImGui::Button("Take screenshot", ImVec2{ITEM_WIDTH, 0}))
    {
        saveScreenshot(screenshotNum++);
    }

    ImGui::Checkbox("Write OBJs", &saveOBJs);
    ImGui::SameLine(ITEM_WIDTH, 2 * INDENT);
    if ((saveOBJs && objNum == 0) || ImGui::Button("Write OBJ", ImVec2{ITEM_WIDTH, 0}))
    {
        saveOBJ(MainApp::instance->mesh, MainApp::instance->geom, MainApp::instance->geomOrig, objNum++);
    }
    ImGui::Checkbox("Log Energy", &MainApp::instance->logEnergy);
    ImGui::SameLine(ITEM_WIDTH, 2 * INDENT);
    ImGui::Checkbox("Topology Change", &MainApp::instance->topoChange);
    ImGui::Checkbox("Show area ratios", &areaRatios);

    const GradientMethod methods[] = {GradientMethod::HsProjectedIterative,
                                      GradientMethod::HsProjected,
                                      GradientMethod::HsExactProjected,
                                      GradientMethod::H1Projected,
                                      GradientMethod::L2Unconstrained,
                                      GradientMethod::L2Projected,
                                      GradientMethod::AQP,
                                      GradientMethod::H1_LBFGS,
                                      GradientMethod::BQN_LBFGS,
                                      GradientMethod::H2Projected};

    selectFromDropdown("Method", methods, IM_ARRAYSIZE(methods), MainApp::instance->methodChoice);

    ImGui::Checkbox("Dynamic remeshing", &remesh);

    const remeshing::RemeshingMode rModes[] = {remeshing::RemeshingMode::FlipOnly,
                                               remeshing::RemeshingMode::SmoothOnly,
                                               remeshing::RemeshingMode::SmoothAndFlip,
                                               remeshing::RemeshingMode::SmoothFlipAndCollapse};

    const remeshing::SmoothingMode sModes[] = {remeshing::SmoothingMode::Laplacian,
                                               remeshing::SmoothingMode::Circumcenter};

    const remeshing::FlippingMode fModes[] = {remeshing::FlippingMode::Delaunay,
                                              remeshing::FlippingMode::Degree};

    selectFromDropdown("Remeshing mode", rModes, IM_ARRAYSIZE(rModes), MainApp::instance->remesher.remeshingMode);
    selectFromDropdown("Smoothing mode", sModes, IM_ARRAYSIZE(sModes), MainApp::instance->remesher.smoothingMode);
    selectFromDropdown("Flipping mode", fModes, IM_ARRAYSIZE(fModes), MainApp::instance->remesher.flippingMode);

    ImGui::Checkbox("Curvature adaptive remeshing", &MainApp::instance->remesher.curvatureAdaptive);

    biorsurfaces::MainApp::instance->HandlePicking();

    ImGui::InputInt("Iteration limit", &MainApp::instance->stepLimit);
    ImGui::InputInt("Real time limit (ms)", &MainApp::instance->realTimeLimit);
    ImGui::InputInt("Topology Change steps", &MainApp::instance->topoChangeSteps);

    if (uiNormalizeView != MainApp::instance->normalizeView)
    {
        biorsurfaces::MainApp::instance->normalizeView = uiNormalizeView;
        biorsurfaces::MainApp::instance->updateMeshPositions();
    }
    ImGui::PopItemWidth();

    ImGui::Checkbox("Limit step size", &limitStep);
    ImGui::SliderFloat("Max step size", &maxStep, 0.001, 0.1);
    ImGui::InputFloat("Weight Helfrich", &regWW); // set a float variable

    ImGui::Text("Model Selection");
    ImGui::BeginGroup();
    ImGui::Indent(INDENT);

    const char* modelItems[] = { "BC Model", "ADE Model" };
    int currentModelIndex = static_cast<int>(currentModel);
    if (ImGui::Combo("Model", &currentModelIndex, modelItems, IM_ARRAYSIZE(modelItems)))
    {
        currentModel = static_cast<ModelType>(currentModelIndex);
        // Reset the model-specific parameters here if needed
        if (currentModel == ModelType::BC) {
            ADE = false;
            regCons = 1.05;
            regalpha = 400;
        } else {
            ADE = true;
            regalpha = 2/ geometrycentral::PI; 
            regCons = 1.05;
        }
    }
    ImGui::EndGroup();

    if (ADE)
    {

      ImGui::InputFloat("alpha or K_ADE", &regalpha); // set a float variable 
      ImGui::InputFloat("Change step of nu and dela0", &changeStep); // set a float variable
      biorsurfaces::MainApp::instance->flow->maxStepSize = limitStep ? maxStep : -1.;
      std::string float_str_dela0 = std::to_string(regCons);
      std::string float_str_step = std::to_string(changeStep);

      ImGui::Text("dela0: %s ", float_str_dela0.c_str());

      ImGui::BeginGroup();
      ImGui::Indent(INDENT);
      std::stringstream ss1, ss2, ss3, ss4;

    // Format the string with a float value
      ss1 << "Add " << std::fixed << std::setprecision(3) << changeStep << " dela0";
      ss2 << "Subtract" << std::fixed << std::setprecision(3) << changeStep << " dela0";
      ss3 << "Add " << std::fixed << std::setprecision(3) << changeStep << " RV";
      ss4 << "Subtract " << std::fixed << std::setprecision(3) << changeStep << " RV";

    // Get the formatted string
      std::string formatted_str1 = ss1.str();
      std::string formatted_str2 = ss2.str();
      std::string formatted_str3 = ss3.str();
      std::string formatted_str4 = ss4.str();

      if (ImGui::Button(formatted_str1.c_str(), ImVec2{ITEM_WIDTH, 0}) || run)
      {
          regCons += changeStep;
          std::cout << "New dela0 is: " << regCons << std::endl;
          biorsurfaces::MainApp::instance->flow->energies[1]->ResetTargetsHelfrich(regCons, regWW, regalpha);
      }

      ImGui::SameLine(ITEM_WIDTH, 2 * INDENT);
      if (ImGui::Button(formatted_str2.c_str(), ImVec2{ITEM_WIDTH, 0}) || run)
      {
          regCons -= changeStep;
          std::cout << "New dela0 is: " << regCons << std::endl;
          biorsurfaces::MainApp::instance->flow->energies[1]->ResetTargetsHelfrich(regCons, regWW,regalpha);
      }
      ImGui::EndGroup();

      std::string float_str_rv = std::to_string(regConsArea);
      ImGui::Text("Reduced Volume: %s", float_str_rv.c_str());

      ImGui::BeginGroup();
      ImGui::Indent(INDENT);

      if (ImGui::Button(formatted_str3.c_str(), ImVec2{ITEM_WIDTH, 0}) || run)
      {
          regConsArea += changeStep;
          std::cout << "New Reduced Volume is: " << regConsArea << std::endl;
          biorsurfaces::MainApp::instance->flow->schurConstraints[1].constraint->setTargetValue(regConsArea);
          biorsurfaces::MainApp::instance->flow->ResetAllConstraints();
      }

      ImGui::SameLine(ITEM_WIDTH, 2 * INDENT);
      if (ImGui::Button(formatted_str4.c_str(), ImVec2{ITEM_WIDTH, 0}) || run)
      {
          regConsArea -= changeStep;
          std::cout << "New Reduced Volume is: " << regConsArea << std::endl;
          biorsurfaces::MainApp::instance->flow->schurConstraints[1].constraint->setTargetValue(regConsArea);
          biorsurfaces::MainApp::instance->flow->ResetAllConstraints();
      }
      ImGui::EndGroup();
    } else {

    ImGui::InputFloat("Weight of BC Term", &regalpha); // set a float variable 
    ImGui::InputFloat("Change step of nu and dela", &changeStep); // set a float variable
    //Important note: 2*regCons is due to the fact that the non-dimensional delta0 is 2 times larger than what is 
    //usually used in the literature for bilayer coupled modeling of biomembranes.
    biorsurfaces::MainApp::instance->flow->energies[1]->ResetTargetsHelfrich(2*regCons, regWW,regalpha/4.0);
    biorsurfaces::MainApp::instance->flow->maxStepSize = limitStep ? maxStep : -1.;
    std::string float_str_rmc = std::to_string(regCons);
    std::string float_str_step = std::to_string(changeStep);

    ImGui::Text("Dela: %s ", float_str_rmc.c_str());

    ImGui::BeginGroup();
    ImGui::Indent(INDENT);
    std::stringstream ss1, ss2, ss3, ss4;

    // Format the string with a float value
    ss1 << "Add " << std::fixed << std::setprecision(3) << changeStep << " Dela";
    ss2 << "Subtract" << std::fixed << std::setprecision(3) << changeStep << " Dela";
    ss3 << "Add " << std::fixed << std::setprecision(3) << changeStep << " nu";
    ss4 << "Subtract " << std::fixed << std::setprecision(3) << changeStep << " nu";

    // Get the formatted string
    std::string formatted_str1 = ss1.str();
    std::string formatted_str2 = ss2.str();
    std::string formatted_str3 = ss3.str();
    std::string formatted_str4 = ss4.str();

    if (ImGui::Button(formatted_str1.c_str(), ImVec2{ITEM_WIDTH, 0}) || run)
    {
        regCons += changeStep;
        std::cout << "New Reduced Mean Curvature is: " << regCons << std::endl;
        biorsurfaces::MainApp::instance->flow->energies[1]->ResetTargetsHelfrich(2*regCons, regWW, regalpha/4.0);
    }

    ImGui::SameLine(ITEM_WIDTH, 2 * INDENT);
    if (ImGui::Button(formatted_str2.c_str(), ImVec2{ITEM_WIDTH, 0}) || run)
    {
        regCons -= changeStep;
        std::cout << "New Reduced Mean Curvature is: " << regCons << std::endl;
        biorsurfaces::MainApp::instance->flow->energies[1]->ResetTargetsHelfrich(2*regCons, regWW, regalpha/4.0);
    }
    ImGui::EndGroup();
    std::string float_str_rv = std::to_string(regConsArea);
    ImGui::Text("Reduced Volume: %s", float_str_rv.c_str());
    ImGui::BeginGroup();
    ImGui::Indent(INDENT);

    if (ImGui::Button(formatted_str3.c_str(), ImVec2{ITEM_WIDTH, 0}) || run)
    {
        regConsArea += changeStep;
        std::cout << "New Reduced Volume is: " << regConsArea << std::endl;
        biorsurfaces::MainApp::instance->flow->schurConstraints[1].constraint->setTargetValue(regConsArea);
        biorsurfaces::MainApp::instance->flow->ResetAllConstraints();
    }

    ImGui::SameLine(ITEM_WIDTH, 2 * INDENT);
    if (ImGui::Button(formatted_str4.c_str(), ImVec2{ITEM_WIDTH, 0}) || run)
    {
        regConsArea -= changeStep;
        std::cout << "New Reduced Volume is: " << regConsArea << std::endl;
        biorsurfaces::MainApp::instance->flow->schurConstraints[1].constraint->setTargetValue(regConsArea);
        biorsurfaces::MainApp::instance->flow->ResetAllConstraints();
    }
    ImGui::EndGroup();
    }

    ImGui::Text("Single Run");

    if (ImGui::Button("Take 1 step", ImVec2{ITEM_WIDTH, 0}) || run)
    {
        MainApp::instance->TakeOptimizationStep(remesh, areaRatios);

        if (takeScreenshots)
        {
            saveScreenshot(screenshotNum++);
        }
        if (saveOBJs)
        {
            saveOBJ(MainApp::instance->mesh, MainApp::instance->geom, MainApp::instance->geomOrig, objNum++);
        }
        if ((MainApp::instance->stepLimit > 0 && MainApp::instance->numSteps >= MainApp::instance->stepLimit) ||
            (MainApp::instance->realTimeLimit > 0 && MainApp::instance->timeSpentSoFar >= MainApp::instance->realTimeLimit))
        {
            run = false;
            if (MainApp::instance->exitWhenDone)
            {
                std::exit(0);
            }
        }
    }


    // ImGui::Text("Remeshing tests");

    // ImGui::BeginGroup();
    // ImGui::Indent(INDENT);

    // if (ImGui::Button("Remesh"))
    // {
    //     MainApp::instance->remesher.Remesh(5, true);
    //     MainApp::instance->mesh->compress();
    //     MainApp::instance->reregisterMesh();
    // }
    // ImGui::EndGroup();

    ImGui::Text("Plotting Gradients");
    ImGui::BeginGroup();
    ImGui::Indent(INDENT);

    if (ImGui::Button("Plot Gradient"))
    {
        MainApp::instance->PlotGradients();
    }
    ImGui::EndGroup();
}

struct MeshAndEnergy
{
    biorsurfaces::TPEKernel *kernel;
    polyscope::SurfaceMesh *psMesh;
    biorsurfaces::MeshPtr mesh;
    biorsurfaces::GeomPtr geom;
    biorsurfaces::UVDataPtr uvs;
    std::string meshName;
};

MeshAndEnergy initTPEOnMesh(std::string meshFile, double alpha, double beta)
{
    using namespace biorsurfaces;
    std::cout << "Initializing tangent-point energy with (" << alpha << ", " << beta << ")" << std::endl;

    MeshUPtr u_mesh;
    std::unique_ptr<VertexPositionGeometry> u_geometry;
    std::unique_ptr<CornerData<Vector2>> uvs;

    // Load mesh
    std::tie(u_mesh, u_geometry, uvs) = readParameterizedMesh(meshFile);
    std::string mesh_name = polyscope::guessNiceNameFromPath(meshFile);

    std::cout << "Read " << uvs->size() << " UV coordinates" << std::endl;
    bool hasUVs = false;

    for (GCVertex v : u_mesh->vertices())
    {
        for (surface::Corner c : v.adjacentCorners())
        {
            Vector2 uv = (*uvs)[c];
            if (uv.x > 0 || uv.y > 0)
            {
                hasUVs = true;
            }
        }
    }

    if (hasUVs)
    {
        std::cout << "Mesh has nonzero UVs; using as flags for attractors" << std::endl;
    }
    else
    {
        std::cout << "Mesh has no UVs or all UVs are 0; not using as flags" << std::endl;
    }

    // Register the mesh with polyscope
    polyscope::SurfaceMesh *psMesh = polyscope::registerSurfaceMesh(mesh_name,
                                                                    u_geometry->inputVertexPositions, u_mesh->getFaceVertexList(),
                                                                    polyscopePermutations(*u_mesh));

    psMesh->setSurfaceColor(glm::vec3(255 / 255., 49 / 255., 49 / 255.));
    psMesh->setEdgeColor(glm::vec3(0.0 / 255., 0.0/ 255., 0.0 / 255.));
    psMesh->setEdgeWidth(1.5);
    psMesh->setSmoothShade(true);
    MeshPtr meshShared = std::move(u_mesh);
    GeomPtr geomShared = std::move(u_geometry);
    UVDataPtr uvShared = std::move(uvs);

    geomShared->requireFaceNormals();
    geomShared->requireFaceAreas();
    geomShared->requireVertexNormals();
    geomShared->requireVertexDualAreas();
    geomShared->requireVertexGaussianCurvatures();

    TPEKernel *tpe = new biorsurfaces::TPEKernel(meshShared, geomShared, alpha, beta);

    std::cout << "Initial mesh area = " << totalArea(geomShared, meshShared) << std::endl;
    std::cout << "Initial mesh volume = " << totalVolume(geomShared, meshShared) << std::endl;

    return MeshAndEnergy{tpe, psMesh, meshShared, geomShared, (hasUVs) ? uvShared : 0, mesh_name};
}

enum class EnergyOverride
{
    TangentPoint,
    Coulomb,
    Helfrich
};

biorsurfaces::SurfaceFlow *setUpFlow(MeshAndEnergy &m, double theta, biorsurfaces::scene::SceneData &scene, EnergyOverride eo)
{
    using namespace biorsurfaces;

    SurfaceEnergy *energy;

    if (eo == EnergyOverride::Coulomb)
    {
        std::cout << "Using Coulomb energy in place of tangent-point energy" << std::endl;
        energy = new CoulombEnergy(m.kernel, theta);
    }
    else if (eo == EnergyOverride::Helfrich)
    {
        std::cout << "Using Helfrich energy in place of tangent-point energy" << std::endl;
        energy = new HelfrichEnergy(m.mesh, m.geom);
    }
    else
    {
        if (theta <= 0)
        {
            std::cout << "Theta was zero (or negative); using exact all-pairs energy." << std::endl;
            energy = new TPEnergyAllPairs(m.kernel->mesh, m.kernel->geom, m.kernel->alpha, m.kernel->beta);
            ;
        }
        else
        {
            std::cout << "Using Barnes-Hut energy with theta = " << theta << "." << std::endl;
            TPEnergyBarnesHut0 *bh = new TPEnergyBarnesHut0(m.kernel->mesh, m.kernel->geom, m.kernel->alpha, m.kernel->beta, theta);
            energy = bh;
        }
    }

    SurfaceFlow *flow = new SurfaceFlow(energy);
    bool kernelRemoved = false;
    flow->allowBarycenterShift = scene.allowBarycenterShift;
    // Set these up here, so that we can aggregate all vertex pins into the same constraint
    Constraints::VertexPinConstraint *pinC = 0;
    Constraints::VertexNormalConstraint *normC = 0;

    std::vector<Vector3> pinLocations;

    for (scene::ConstraintData &data : scene.constraints)
    {
        switch (data.type)
        {
        case scene::ConstraintType::Barycenter:
            kernelRemoved = true;
            flow->addSimpleConstraint<Constraints::BarycenterConstraint3X>(m.mesh, m.geom);
            break;
        case scene::ConstraintType::TotalArea:
            flow->addSchurConstraint<Constraints::TotalAreaConstraint>(m.mesh, m.geom, data.targetMultiplier, data.numIterations, data.targetAddition);
            break;
        case scene::ConstraintType::TotalVolume:
            flow->addSchurConstraint<Constraints::TotalVolumeConstraint>(m.mesh, m.geom, data.targetMultiplier, data.numIterations, data.targetAddition);
            break;
        case scene::ConstraintType::ReducedVolume:
            flow->addSchurConstraint<Constraints::ReducedVolume>(m.mesh, m.geom, data.targetMultiplier, data.numIterations, data.targetAddition);
            break;
        case scene::ConstraintType::ReducedMC:
            flow->addSchurConstraint<Constraints::ReducedMC>(m.mesh, m.geom, data.targetMultiplier, data.numIterations, data.targetAddition);
            break;
        case scene::ConstraintType::BoundaryPins:
        {
            if (!pinC)
            {
                pinC = flow->addSimpleConstraint<Constraints::VertexPinConstraint>(m.mesh, m.geom);
            }
            // Manually add all of the boundary vertex indices as pins
            std::vector<size_t> boundaryInds;
            VertexIndices inds = m.mesh->getVertexIndices();
            for (GCVertex v : m.mesh->vertices())
            {
                if (v.isBoundary())
                {
                    boundaryInds.push_back(inds[v]);

                    Vector3 pos = m.geom->inputVertexPositions[v];
                    pinLocations.push_back(pos);
                }
            }
            pinC->pinVertices(m.mesh, m.geom, boundaryInds);
            kernelRemoved = true;
        }
        break;

        case scene::ConstraintType::VertexPins:
        {
            if (!pinC)
            {
                pinC = flow->addSimpleConstraint<Constraints::VertexPinConstraint>(m.mesh, m.geom);
            }
            // Add the specified vertices as pins
            pinC->pinVertices(m.mesh, m.geom, scene.vertexPins);
            for (VertexPinData &pinData : scene.vertexPins)
            {
                Vector3 pos = m.geom->inputVertexPositions[pinData.vertID];
                pinLocations.push_back(pos);
            }
            // Clear the data vector so that we don't add anything twice
            scene.vertexPins.clear();
            kernelRemoved = true;
        }
        break;

        case scene::ConstraintType::BoundaryNormals:
        {
            if (!normC)
            {
                normC = flow->addSimpleConstraint<Constraints::VertexNormalConstraint>(m.mesh, m.geom);
            }
            // Manually add all of the boundary vertex indices as pins
            std::vector<size_t> boundaryInds;
            VertexIndices inds = m.mesh->getVertexIndices();
            for (GCVertex v : m.mesh->vertices())
            {
                if (v.isBoundary())
                {
                    boundaryInds.push_back(inds[v]);
                }
            }
            normC->pinVertices(m.mesh, m.geom, boundaryInds);
        }

        case scene::ConstraintType::VertexNormals:
        {
            if (!normC)
            {
                normC = flow->addSimpleConstraint<Constraints::VertexNormalConstraint>(m.mesh, m.geom);
            }
            // Add the specified vertices as pins
            normC->pinVertices(m.mesh, m.geom, scene.vertexNormals);
            // Clear the data vector so that we don't add anything twice
            scene.vertexNormals.clear();
        }
        break;

        default:
            std::cout << "  * Skipping unrecognized constraint type" << std::endl;
            break;
        }
    }

    if (!kernelRemoved)
    {
        // std::cout << "Auto-adding barycenter constraint to eliminate constant kernel of Laplacian" << std::endl;
        // flow->addSimpleConstraint<Constraints::BarycenterConstraint3X>(m.mesh, m.geom);
    }

    if (pinLocations.size() > 0)
    {
        polyscope::registerPointCloud("pinned vertices", pinLocations);
    }

    return flow;
}

biorsurfaces::scene::SceneData defaultScene(std::string meshName)
{
    using namespace biorsurfaces;
    using namespace biorsurfaces::scene;
    SceneData data;
    data.meshName = meshName;
    data.alpha = 6;
    data.beta = 12;
    data.constraints = std::vector<ConstraintData>({ConstraintData{scene::ConstraintType::Barycenter, 1, 0, 0},
                                                    ConstraintData{scene::ConstraintType::TotalArea, 1, 0, 0},
                                                    ConstraintData{scene::ConstraintType::TotalVolume, 1, 0, 0}});
    return data;
}

int main(int argc, char **argv)
{
    using namespace biorsurfaces;

    // Configure the argument parser
    args::ArgumentParser parser("Repulsive biomembranes: A tool for Discrete Differential Analysis of Biomembranes");
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");
    args::ValueFlag<double> thetaFlag(parser, "Theta", "Theta value for Barnes-Hut approximation; 0 means exact.", args::Matcher{'t', "theta"});
    args::ValueFlag<std::string> mult_alg_Flag(parser, "mult_alg", "Algorithm for the near field matrix-vector product. Possible values are \"Hybrid\" (default) and \"MKL_CSR\" (maybe more robust).", {"mult_alg"});
    args::ValueFlagList<std::string> obstacleFiles(parser, "obstacles", "Obstacles to add", {'o'});
    args::Flag autologFlag(parser, "autolog", "Automatically start the flow, log performance, and exit when done.", {"autolog"});
    args::Flag coulombFlag(parser, "coulomb", "Use a coulomb energy instead of the tangent-point energy.", {"coulomb"});
    args::ValueFlag<int> threadFlag(parser, "threads", "How many threads to use in parallel.", {"threads"});
    args::ValueFlag<double> areaFlag(parser, "ra", "Reduced area", args::Matcher{"ra"});
    args::ValueFlag<double> volFlag(parser, "rv", "Reduced volume", args::Matcher{"rv"});

    polyscope::options::programName = "Repulsive BioSurfaces- Meghdad Razizadeh(2023) -SJCRH";
    polyscope::options::groundPlaneEnabled = false;

    std::cout << "Using Eigen version " << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION << std::endl;

    MKLVersion Version;
    mkl_get_version(&Version);

    std::cout << "Using MKL version " << Version.MajorVersion << "." << Version.MinorVersion << "." << Version.UpdateVersion << std::endl;

    // Parse args
    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help)
    {
        std::cout << parser;
        return 0;
    }
    catch (args::ParseError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    // Make sure a mesh name was given
    if (!inputFilename)
    {
        std::cerr << "Please specify a mesh file as argument" << std::endl;
        return EXIT_FAILURE;
    }

    MainApp::defaultNumThreads = omp_get_max_threads() / 2 + 2;

    if (threadFlag)
    {
        int nThreads = args::get(threadFlag);
        std::cout << "Using " << nThreads << " threads as specified." << std::endl;
        omp_set_num_threads(nThreads);
        MainApp::specifiedNumThreads = nThreads;
    }
    else
    {
        omp_set_num_threads(MainApp::defaultNumThreads);
        MainApp::specifiedNumThreads = MainApp::defaultNumThreads;
        std::cout << "Defaulting to " << MainApp::defaultNumThreads << " threads." << std::endl;
    }

    if (mult_alg_Flag)
    {
        std::string s = args::get(mult_alg_Flag);
        if (s.compare("MKL_CSR") == 0)
        {
            BCTDefaultSettings.mult_alg = NearFieldMultiplicationAlgorithm::MKL_CSR;
            std::cout << "Using \"MKL_CSR\" for near field matrix-vector product." << std::endl;
        }
        else if (s.compare("Hybrid") == 0)
        {
            BCTDefaultSettings.mult_alg = NearFieldMultiplicationAlgorithm::Hybrid;
            std::cout << "Using \"Hybrid\" for near field matrix-vector product." << std::endl;
        }
        else
        {
            BCTDefaultSettings.mult_alg = NearFieldMultiplicationAlgorithm::Hybrid;
            std::cout << "Unknown method \"" + s + "\". Using default value \"Hybrid\" for near field matrix-vector product." << std::endl;
        }
    }
    else
    {
        BCTDefaultSettings.mult_alg = NearFieldMultiplicationAlgorithm::Hybrid;
        std::cout << "Using default value \"Hybrid\" for near field matrix-vector product." << std::endl;
    }

    double theta = 0.3;
    if (!thetaFlag)
    {
        std::cout << "Barnes-Hut theta value not specified; defaulting to theta = " << theta << std::endl;
    }
    else
    {
        theta = args::get(thetaFlag);
    }

    // Initialize polyscope
    polyscope::init();
    // Set the callback function
    polyscope::state::userCallback = customCallback;

    // Parse the input file, either as a scene file or as a mesh
    std::string inFile = args::get(inputFilename);
    scene::SceneData data;

    if (endsWith(inFile, ".txt") || endsWith(inFile, ".scene"))
    {
        std::cout << "Parsing " << inFile << " as scene file." << std::endl;
        data = scene::parseScene(inFile);
    }

    else if (endsWith(inFile, ".obj"))
    {
        std::cout << "Parsing " << inFile << " as OBJ mesh file." << std::endl;
        data = defaultScene(inFile);
    }

    else
    {
        throw std::runtime_error("Unknown file extension for " + inFile + ".");
    }

    bool useCoulomb = false;
    if (coulombFlag)
    {
        useCoulomb = true;
        std::cout << "Using Coulomb energy. (Note: Not expected to work well.)" << std::endl;
    }

    MeshAndEnergy m = initTPEOnMesh(data.meshName, data.alpha, data.beta);

    EnergyOverride eo = EnergyOverride::TangentPoint;
    if (useCoulomb)
    {
        eo = EnergyOverride::Coulomb;
    }
    else if (data.defaultMethod == GradientMethod::Helfrich)
    {
        eo = EnergyOverride::Helfrich;
    }

    SurfaceFlow *flow = setUpFlow(m, theta, data, eo);
    flow->disableNearField = data.disableNearField;

    MainApp::instance = new MainApp(m.mesh, m.geom, flow, m.psMesh, m.meshName);
    MainApp::instance->bh_theta = theta;
    MainApp::instance->kernel = m.kernel;
    MainApp::instance->stepLimit = data.iterationLimit;
    MainApp::instance->realTimeLimit = data.realTimeLimit;
    MainApp::instance->methodChoice = data.defaultMethod;
    MainApp::instance->sceneData = data;
    MainApp::instance->uvs = m.uvs;

    if (autologFlag)
    {
        std::cout << "Autolog flag was used; starting flow automatically." << std::endl;
        MainApp::instance->exitWhenDone = true;
        MainApp::instance->logEnergy = true;
        MainApp::instance->topoChange = true;
        MainApp::instance->topoChangeSteps = 5;
        run = true;
        std::ofstream outfile;
        outfile.open(data.performanceLogFile, std::ios_base::out);
        outfile.close();
    }

    for (scene::PotentialData &p : data.potentials)
    {
        MainApp::instance->AddPotential(p.type, p.weight, p.targetValue);
    }
    for (scene::ObstacleData &obs : data.obstacles)
    {
        MainApp::instance->AddObstacle(obs.obstacleName, obs.weight, obs.recenter, obs.asPointCloud);
    }
    for (scene::ImplicitBarrierData &barrierData : data.implicitBarriers)
    {
        MainApp::instance->AddImplicitBarrier(barrierData);
    }

    if (data.autoComputeVolumeTarget)
    {
        double targetVol = MainApp::instance->totalObstacleVolume * data.autoVolumeTargetRatio;
        std::cout << "Retargeting volume constraint to value " << targetVol << " (" << data.autoVolumeTargetRatio << "x obstacle volume)" << std::endl;
        MainApp::instance->flow->retargetSchurConstraintOfType<Constraints::TotalVolumeConstraint>(targetVol);
    }

    MainApp::instance->updateMeshPositions();
    MainApp::instance->remesher.curvatureAdaptive = false;
    float regWW = 5000;
    float changeStep = 0.025;
    MainApp::instance->flow->energies[1]->ResetTargetsHelfrich(1.1, regWW,1.0);
    bool remesh = true;
    bool areaRatios = false;
    uint objNum = 0;
    uint screenshotNum = 0;
    std::ofstream logFile;
    logFile.open("log_Energy.txt", std::ios::trunc);
    logFile << "HelfrichEnergy"
            << " "
            << "TPEEnergy"
            << " "
            << std::endl;
    if (!logFile.is_open())
    {
        std::cerr << "Error opening file" << std::endl;
        return 1;
    }
    
    polyscope::show();
    return EXIT_SUCCESS;
}
