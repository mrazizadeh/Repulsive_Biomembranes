#include "sobolev/h1.h"
#include "helpers.h"
#include "sobolev/all_constraints.h"

namespace biorsurfaces
{
    namespace H1
    {
        void getTriplets(std::vector<Triplet> &triplets, const MeshPtr &mesh, const GeomPtr &geom, double epsilon, bool premultiplyMass)
        {
            geom->requireEdgeCotanWeights();
            VertexIndices indices = mesh->getVertexIndices();
            for (GCVertex v : mesh->vertices())
            {
                double rowSum = 0;
                double area = geom->vertexDualAreas[v];
                double mass = premultiplyMass ? geom->vertexDualAreas[v] : 1;

                for (GCEdge e : v.adjacentEdges())
                {
                    GCVertex opp = getOppositeVertex(e, v);
                    double wt = geom->edgeCotanWeight(e);
                    rowSum += wt;
                    triplets.push_back(Triplet(indices[v], indices[opp], -wt / mass));
                }
                triplets.push_back(Triplet(indices[v], indices[v], rowSum / mass + epsilon * geom->vertexDualAreas[v]));
            }
            geom->unrequireEdgeCotanWeights();
        }

        void ProjectGradient(Eigen::MatrixXd &gradient, Eigen::MatrixXd &dest, MeshPtr &mesh, GeomPtr &geom, bool useMass)
        {
            // Assemble the metric matrix
            std::vector<Triplet> triplets, triplets3x;
            H1::getTriplets(triplets, mesh, geom, false);
            // Reduplicate the entries 3x along diagonals
            MatrixUtils::TripleTriplets(triplets, triplets3x);
            Constraints::BarycenterConstraint3X bconstraint(mesh, geom);
            Constraints::addTripletsToSymmetric(bconstraint, triplets3x, mesh, geom, 3 * mesh->nVertices());
            Eigen::SparseMatrix<double> metric(3 * mesh->nVertices() + 3, 3 * mesh->nVertices() + 3);
            metric.setFromTriplets(triplets3x.begin(), triplets3x.end());

            dest = gradient;
            if (useMass)
            {
                MultiplyByMass(dest, mesh, geom);
            }

            // Flatten the gradient into a single column
            Eigen::VectorXd gradientCol;
            gradientCol.setZero(3 * mesh->nVertices() + 3);
            MatrixUtils::MatrixIntoColumn(dest, gradientCol);

            // Invert the metric, and write it into the destination
            MatrixUtils::SolveSparseSystem(metric, gradientCol, gradientCol);
            MatrixUtils::ColumnIntoMatrix(gradientCol, dest);
        }

        void ProjectConstraints(MeshPtr &mesh, GeomPtr &geom, std::vector<Constraints::SimpleProjectorConstraint *> simpleConstraints,
                                std::vector<ConstraintPack> &newtonConstraints, SparseFactorization &factoredL, int newtonIterations)
        {
            size_t nRows = factoredL.nRows;
            Eigen::VectorXd errors(nRows);

            int nIters = 0;

            while (nIters < newtonIterations)
            {
                // Fill right-hand side with error values
                errors.setZero();

                size_t curRow = 3 * mesh->nVertices();
                for (Constraints::SimpleProjectorConstraint *spc : simpleConstraints)
                {
                    spc->addErrorValues(errors, mesh, geom, curRow);
                    curRow += spc->nRows();
                }
                for (ConstraintPack &c : newtonConstraints)
                {
                    c.constraint->addErrorValues(errors, mesh, geom, curRow);
                    curRow += c.constraint->nRows();
                }

                double maxError = errors.lpNorm<Eigen::Infinity>();
                std::cout << "Error after " << nIters << " iterations = " << maxError << std::endl;
                if (maxError < 1e-3)
                {
                    break;
                }
                nIters++;

                errors = factoredL.Solve(errors);

                // Apply the correction
                for (size_t i = 0; i < mesh->nVertices(); i++)
                {
                    GCVertex v_i = mesh->vertex(i);
                    Vector3 corr_i{errors(3 * i), errors(3 * i + 1), errors(3 * i + 2)};
                    geom->inputVertexPositions[v_i] -= corr_i;
                }
            }
        }
    } // namespace H1
} // namespace biorsurfaces