#include "sobolev/constraints/vertex_normal.h"
#include "helpers.h"
#include "surface_derivatives.h"

namespace biorsurfaces
{
    namespace Constraints
    {
        VertexNormalConstraint::VertexNormalConstraint(const MeshPtr &mesh, const GeomPtr &geom) {}

        void VertexNormalConstraint::pinVertices(const MeshPtr &mesh, const GeomPtr &geom, std::vector<size_t> &indices_)
        {
            for (size_t i = 0; i < indices_.size(); i++)
            {
                indices.push_back(indices_[i]);
                initNormals.push_back(vertexAreaNormal(geom, mesh->vertex(indices_[i])));
            }
        }

        void VertexNormalConstraint::ResetFunction(const MeshPtr &mesh, const GeomPtr &geom)
        {
            for (size_t i = 0; i < indices.size(); i++)
            {
                initNormals[i] = vertexAreaNormal(geom, mesh->vertex(indices[i]));
            }
        }

        void VertexNormalConstraint::addTriplets(std::vector<Triplet> &triplets, const MeshPtr &mesh, const GeomPtr &geom, int baseRow)
        {
            VertexIndices allInds = mesh->getVertexIndices();
            // All we do is put a 1 in the index for all pinned vertices
            for (size_t i = 0; i < indices.size(); i++)
            {
                size_t thisRow = baseRow + 3 * i;
                GCVertex v_i = mesh->vertex(indices[i]);

                // Need to add Jacobians for all neighbors
                for (GCVertex neighbor : v_i.adjacentVertices())
                {
                    Jacobian J_neighbor = SurfaceDerivs::vertexNormalWrtVertex(geom, v_i, neighbor);
                    size_t neighborCol = 3 * allInds[neighbor];
                    J_neighbor.AddTransposeTriplets(triplets, thisRow, neighborCol);
                }

                // Also add Jacobian for center vertex
                Jacobian J_self = SurfaceDerivs::vertexNormalWrtVertex(geom, v_i, v_i);
                size_t thisCol = 3 * allInds[v_i];
                J_self.AddTransposeTriplets(triplets, thisRow, thisCol);
            }
        }

        void VertexNormalConstraint::addEntries(Eigen::MatrixXd &M, const MeshPtr &mesh, const GeomPtr &geom, int baseRow)
        {
            VertexIndices allInds = mesh->getVertexIndices();
            // All we do is put a 1 in the index for all pinned vertices
            for (size_t i = 0; i < indices.size(); i++)
            {
                size_t thisRow = baseRow + 3 * i;
                GCVertex v_i = mesh->vertex(indices[i]);

                // Need to add Jacobians for all neighbors
                for (GCVertex neighbor : v_i.adjacentVertices())
                {
                    Jacobian J_neighbor = SurfaceDerivs::vertexNormalWrtVertex(geom, v_i, neighbor);
                    size_t neighborCol = 3 * allInds[neighbor];
                    J_neighbor.AddTransposeToMatrix(M, thisRow, neighborCol);
                }

                // Also add Jacobian for center vertex
                Jacobian J_self = SurfaceDerivs::vertexNormalWrtVertex(geom, v_i, v_i);
                size_t thisCol = 3 * allInds[v_i];
                J_self.AddTransposeToMatrix(M, thisRow, thisCol);
            }
        }

        void VertexNormalConstraint::addErrorValues(Eigen::VectorXd &V, const MeshPtr &mesh, const GeomPtr &geom, int baseRow)
        {
            for (size_t i = 0; i < indices.size(); i++)
            {
                Vector3 current = vertexAreaNormal(geom, mesh->vertex(indices[i]));
                V(baseRow + 3 * i) = (current.x - initNormals[i].x);
                V(baseRow + 3 * i + 1) = (current.y - initNormals[i].y);
                V(baseRow + 3 * i + 2) = (current.z - initNormals[i].z);
            }
        }

        size_t VertexNormalConstraint::nRows()
        {
            // Every vertex gets 3 rows
            return indices.size() * 3;
        }

        void VertexNormalConstraint::ProjectConstraint(MeshPtr &mesh, GeomPtr &geom)
        {
            return;
        }

    } // namespace Constraints
} // namespace biorsurfaces