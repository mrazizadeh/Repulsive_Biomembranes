#include "deprecated/barnes_hut_tpe_6d.h"
#include "helpers.h"
#include "surface_derivatives.h"
#include "matrix_utils.h"

namespace biorsurfaces
{

    BarnesHutTPEnergy6D::BarnesHutTPEnergy6D(TPEKernel *kernel_, double theta_)
    {
        kernel = kernel_;
        theta = theta_;
        root = 0;
        disableNearField = false;
    }

    BarnesHutTPEnergy6D::~BarnesHutTPEnergy6D()
    {
        if (root)
        {
            delete root;
        }
    }

    double BarnesHutTPEnergy6D::Value()
    {
        kernel->recomputeBarycenters();

        if (!root)
        {
            std::cerr << "BVH for BarnesHutTPEnergy6D was not initialized. Call Update() first." << std::endl;
            std::exit(1);
        }

        double sum = 0;
        #pragma omp parallel for reduction(+ : sum) shared(root, energyPerFace)
        for (size_t i = 0; i < kernel->mesh->nFaces(); i++)
        {
            GCFace f = kernel->mesh->face(i);
            double e_f = computeEnergyOfFace(f, root);
            energyPerFace[f] = e_f;
            sum += e_f;
        }
        return sum;
    }

    double BarnesHutTPEnergy6D::computeEnergyOfFace(GCFace face, BVHNode6D *bvhRoot)
    {
        Vector3 bcenter = kernel->faceBarycenters[face];

        if (bvhRoot->nodeType == BVHNodeType::Empty)
        {
            // Empty cluster = no value
            return 0;
        }
        else if (bvhRoot->nodeType == BVHNodeType::Leaf)
        {
            // If we disable the near field, then these get skipped
            if (false)
            {
                return 0;
            }
            else
            {
                // Compute the energy exactly for the one face in the cluster
                GCFace f2 = bvhRoot->getSingleFace(kernel->mesh);
                return kernel->tpe_pair(face, f2);
            }
        }
        else if (bvhRoot->isAdmissibleFrom(bcenter, theta))
        {
            // Use the cluster approximation
            MassNormalPoint mnp = bvhRoot->GetMassNormalPoint();
            return kernel->tpe_pair(face, mnp);
        }
        else
        {
            // Recursively compute it on all children
            double sum = 0;
            for (BVHNode6D *child : bvhRoot->children)
            {
                sum += computeEnergyOfFace(face, child);
            }
            return sum;
        }
    }

    void BarnesHutTPEnergy6D::Differential(Eigen::MatrixXd &output)
    {
        if (!root)
        {
            std::cerr << "BVH for BarnesHutTPEnergy6D was not initialized. Call Update() first." << std::endl;
            std::exit(1);
        }

        kernel->recomputeBarycenters();

        VertexIndices indices = kernel->mesh->getVertexIndices();
        output.setZero();
        Eigen::MatrixXd partialOutput = output;
        #pragma omp parallel firstprivate(partialOutput) shared(root, output)
        {
            #pragma omp for
            for (size_t i = 0; i < kernel->mesh->nFaces(); i++)
            {
                GCFace f = kernel->mesh->face(i);
                accumulateTPEGradient(partialOutput, root, f, indices);
            }

            #pragma omp critical
            {
                output += partialOutput;
            }
        }
    }

    // Add derivatives of all energy terms of the form (f1, _) or (_, f1)
    // with respect to the neighbor vertices of f1.
    void BarnesHutTPEnergy6D::accumulateTPEGradient(Eigen::MatrixXd &gradients, BVHNode6D *node, GCFace face1,
                                                    surface::VertexData<size_t> &indices)
    {
        if (node->nodeType == BVHNodeType::Empty)
        {
            return;
        }
        else if (node->nodeType == BVHNodeType::Leaf)
        {
            // If near field is disabled, then these get skipped
            if (false)
            {
                return;
            }

            // If this is a leaf, then it only has one face in it, so just use it
            GCFace face2 = node->getSingleFace(kernel->mesh);
            // Skip if the faces are the same
            if (face1 == face2)
                return;

            // Differentiate by adjacent vertices to face1
            std::vector<GCVertex> neighborVerts;

            for (GCVertex v : face1.adjacentVertices())
            {
                // Add the forward term (f1, f2)
                Vector3 deriv1 = kernel->tpe_gradient_pair(face1, face2, v);
                size_t r = indices[v];
                MatrixUtils::addToRow(gradients, r, deriv1);

                // Determine if the reverse term (f2, f1) should be added.
                // If v is also adjacent to f2, then it shouldn't be, because
                // (f2, f1) wrt v will be accumulated later.
                bool noOverlap = true;
                for (GCVertex v2 : face2.adjacentVertices())
                {
                    if (v == v2)
                    {
                        noOverlap = false;
                        break;
                    }
                }
                // If v is not adjacent to v2, then we need to add (f2, f1) wrt v now.
                if (noOverlap)
                {
                    Vector3 deriv2 = kernel->tpe_gradient_pair(face2, face1, v);
                    MatrixUtils::addToRow(gradients, r, deriv2);
                }
            }
        }
        else
        {
            Vector3 f1_center = kernel->faceBarycenters[face1];
            if (node->isAdmissibleFrom(f1_center, theta))
            {
                Vector3 normal = node->averageNormal;
                // This cell is far enough away that we can treat it as a single body
                MassNormalPoint mnp2 = node->GetMassNormalPoint();

                // Differentiate both terms for all neighbors
                for (GCVertex v : face1.adjacentVertices())
                {
                    // Derivatives of both forward and reverse terms
                    MatrixUtils::addToRow(gradients, indices[v], kernel->tpe_gradient_pair(face1, mnp2, v));
                    MatrixUtils::addToRow(gradients, indices[v], kernel->tpe_gradient_pair(mnp2, face1, v));
                }
            }
            else
            {
                // Otherwise we continue recursively traversing the tree
                for (size_t i = 0; i < BVH_N_CHILDREN; i++)
                {
                    accumulateTPEGradient(gradients, node->children[i], face1, indices);
                }
            }
        }
    }

    void BarnesHutTPEnergy6D::Update()
    {
        /*
        if (root)
        {
            delete root;
        }

        root = Create6DBVHFromMeshFaces(kernel->mesh, kernel->geom);
        energyPerFace = geometrycentral::surface::FaceData<double>(*kernel->mesh);
        */
       throw std::runtime_error("deprecated");
    }

    MeshPtr BarnesHutTPEnergy6D::GetMesh()
    {
        return kernel->mesh;
    }

    GeomPtr BarnesHutTPEnergy6D::GetGeom()
    {
        return kernel->geom;
    }

    Vector2 BarnesHutTPEnergy6D::GetExponents()
    {
        return Vector2{kernel->alpha, kernel->beta};
    }

    OptimizedClusterTree *BarnesHutTPEnergy6D::GetBVH()
    {
        // return root;
        throw std::runtime_error("deprecated");
        return 0;
    }

    double BarnesHutTPEnergy6D::GetTheta()
    {
        return theta;
    }

} // namespace biorsurfaces
