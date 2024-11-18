#pragma once

#include "rsurface_types.h"
#include "matrix_utils.h"
#include "constraints.h"
#include "optimized_bct.h"
#include "hs_operators.h"
#include "energy/tp_obstacle_barnes_hut_0.h"
#include "sobolev/h1.h"
#include "sobolev/all_constraints.h"
#include "sobolev/sparse_factorization.h"
#include "bct_constructors.h"
#include "metric_term.h"

#include <Eigen/Sparse>

namespace biorsurfaces
{
    using Constraints::ConstraintBase;

    namespace Hs
    {
        struct SchurComplement
        {
            Eigen::MatrixXd C;
            Eigen::MatrixXd M_A;
            Eigen::MatrixXd Ainv_CT;
        };

        Vector3 HatGradientOnTriangle(GCFace face, GCVertex vert, GeomPtr &geom);
        double get_s(double alpha, double beta);

        class HsMetric;

        template <typename Inverse, typename HsPtr>
        void GetSchurComplement(const HsPtr hs, SchurComplement &dest);

        void ProjectViaSchurV(const HsMetric &hs, Eigen::VectorXd &gradient, Eigen::VectorXd &dest);

        /*
        * This class contains everything necessary to compute an Hs-projected
        * gradient direction. A new instance should be used every timestep.
        */
        class HsMetric
        {
        public:
            HsMetric(std::vector<SurfaceEnergy*> energies, SurfaceEnergy *obstacleEnergy_,
                     std::vector<Constraints::SimpleProjectorConstraint *> &spcs,
                     std::vector<ConstraintPack> &schurs);
            ~HsMetric();

            // Build the "high order" fractional Laplacian of order 2s.
            void FillMatrixHigh(Eigen::MatrixXd &M, double s, const MeshPtr &mesh, const GeomPtr &geom) const;
            // Add the regularizing "low order" term.
            void FillMatrixLow(Eigen::MatrixXd &M, double s, const MeshPtr &mesh, const GeomPtr &geom) const;

            // Build the base fractional Laplacian of order s.
            void FillMatrixFracOnly(Eigen::MatrixXd &M, double s, const MeshPtr &mesh, const GeomPtr &geom) const;
            // Build an exact Hs preconditioner with high- and low-order terms.
            Eigen::MatrixXd GetHsMatrixConstrained(bool includeNewton = true) const;

            // Build just the constraint block of the saddle matrix.
            Eigen::SparseMatrix<double> GetConstraintBlock(bool includeNewton = true) const;

            void ProjectGradientExact(const Eigen::MatrixXd &gradient, Eigen::MatrixXd &dest, Eigen::PartialPivLU<Eigen::MatrixXd> &solver) const;

            void ProjectSimpleConstraints();
            void ProjectSimpleConstraintsWithSaddle();

            template <typename Rhs>
            inline Rhs InvertSparseForIterative(const Rhs &gradient) const
            {
                double epsilon = (mesh->nConnectedComponents() > 1) ? 1e-2 : 1e-8;
                Eigen::VectorXd temp = gradient;
                ProjectSparse(temp, temp, epsilon);
                return Rhs(temp);
            }

            template <typename Rhs>
            inline Rhs InvertMetricSchurTemplated(const Rhs &gradient) const
            {
                Eigen::VectorXd temp = gradient;
                ProjectViaSchurV(*this, temp, temp);
                return Rhs(temp);
            }

            template <typename V, typename Dst>
            inline void InvertMetric(const V &gradient, Dst &dest) const
            {
                ProjectSparse(gradient, dest);
            }

            inline void InvertMetricMat(const Eigen::MatrixXd &gradient, Eigen::MatrixXd &dest) const
            {
                ProjectSparseMat(gradient, dest);
            }

            inline double getHsOrder() const
            {
                Vector2 exps = energy->GetExponents();
                return get_s(exps.x, exps.y);
            }

            inline size_t getNumConstraints() const
            {
                return simpleRows + newtonRows;
            }

            inline size_t getNumRows(bool includeNewton = true) const
            {
                if (includeNewton)
                {
                    return 3 * mesh->nVertices() + getNumConstraints();
                }
                else
                {
                    return 3 * mesh->nVertices() + simpleRows;
                }
            }

            inline OptimizedClusterTree *GetBVH() const
            {
                return bvh;
            }

            inline double getBHTheta() const
            {
                return bh_theta;
            }

            inline void SetNewtonConstraints(std::vector<ConstraintPack> &newConstrs)
            {
                newtonConstraints.clear();
                newtonConstraints = newConstrs;
            }

            inline void SetSimpleConstraints(std::vector<Constraints::SimpleProjectorConstraint *> &simples)
            {
                simpleConstraints.clear();
                simpleConstraints = simples;
            }

            inline size_t topLeftNumRows(bool subComponentBarycenters = true) const
            {
                return 3 * mesh->nVertices() + simpleRows;
            }

            mutable MeshPtr mesh;
            mutable GeomPtr geom;
            bool allowBarycenterShift;

            std::vector<Constraints::SimpleProjectorConstraint *> simpleConstraints;
            std::vector<ConstraintPack> newtonConstraints;

            inline void ResetSchurComplement()
            {
                if (schurComplementComputed)
                {
                    schurComplementComputed = false;
                }
            }

            template <typename Inverse>
            inline SchurComplement &Schur() const
            {
                if (!schurComplementComputed)
                {
                    GetSchurComplement<Inverse>(this, schurComplement);
                    schurComplementComputed = true;
                }
                return schurComplement;
            }
            void shiftBarycenterConstraint(Vector3 shift);

            inline std::vector<MetricTerm*>& getMetricTerms() const
            {
                if (metricTerms.size() == 0)
                {
                    // Add the main Hs metric term
                    metricTerms.push_back(new BCTMetricTerm(getBlockClusterTree()));
                    // Can add more metric terms based on extra energies below here
                    for (SurfaceEnergy* extra : extraEnergies)
                    {
                        extra->AddMetricTerm(metricTerms);
                    }
                }
                return metricTerms;
            }

            bool disableNearField = false;

            inline BCTPtr getBlockClusterTree() const
            {
                if (!optBCT)
                {
                    Vector2 exps = energy->GetExponents();
                    // Configuring BCT
                    BCTSettings settings;
                    if (disableNearField)
                    {
                        // This actually turns off both the near and far field of the the low order term. Still somewhat experimental and maybe obsolete.
                        settings.near_lo_modifier = 0.;
                        settings.far_lo_modifier = 0.;
                        std::cout << "    * Low-order near-field interactions in metric BCT are disabled." << std::endl;
                    }
                    // Now this tells the BCT to multiply the metrics by the energy's weight.
                    optBCT = CreateOptimizedBCTFromBVH(bvh, exps.x, exps.y, bh_theta, energy->GetWeight(), settings);


                    if (obstacleEnergy)
                    {
                        std::cout << "    * Building obstacle BCT" << std::endl;
                        OptimizedClusterTree* obstacleBVH = obstacleEnergy->GetBVH();
                        std::cout << "    * Got obstacle BVH " << obstacleBVH << " (" << obstacleBVH->cluster_count << " clusters)" << std::endl;
                        BCTSettings settings;
                        if (disableNearField)
                        {
                            settings.near_lo_modifier = 0.;
                            settings.far_lo_modifier = 0.;
                            std::cout << "    * Low-order near-field interactions in obstacle BCT are disabled." << std::endl;
                        }
                        // Now this tells the BCT to multiply the metrics by the obstacleEnergy's weight.
                        // Should be useful for regularization/penalty scenarios in which one wants to apply extremely large weights.
                        obstacleBCT = std::make_shared<OptimizedBlockClusterTree>(bvh, obstacleBVH, exps.x, exps.y, bh_theta, obstacleEnergy->GetWeight(), settings);
                        std::cout << "    * Built obstacle BCT" << std::endl;
                        optBCT->AddObstacleCorrection(obstacleBCT);
                        std::cout << "    * Added obstacle correction" << std::endl;
                    }

                }
                return optBCT;
            }

        private:
            void addSimpleConstraintEntries(Eigen::MatrixXd &M) const;
            void addSimpleConstraintTriplets(std::vector<Triplet> &triplets) const;
            void initFromEnergy(SurfaceEnergy *energy_);
            void precomputeSizes();

            // Project the gradient into Hs by using the L^{-1} M L^{-1} factorization
            template <typename V, typename Dst>
            void ProjectSparse(const V &gradient, Dst &dest, double epsilon = 1e-10) const;
            // Same as above but with the input/output being matrices
            void ProjectSparseMat(const Eigen::MatrixXd &gradient, Eigen::MatrixXd &dest, double epsilon = 1e-10) const;

            mutable std::vector<MetricTerm*> metricTerms;
            OptimizedClusterTree *bvh;
            double bh_theta;
            double bct_theta;

            size_t simpleRows;
            size_t newtonRows;

            SurfaceEnergy *energy;
            std::vector<SurfaceEnergy*> extraEnergies;

            SurfaceEnergy *obstacleEnergy;
            bool usedDefaultConstraint;

            mutable SparseFactorization factorizedLaplacian;
            mutable BCTPtr optBCT;
            mutable BCTPtr obstacleBCT;
            mutable bool schurComplementComputed;
            mutable SchurComplement schurComplement;
        };

        template <typename V, typename Dst>
        void HsMetric::ProjectSparse(const V &gradientCol, Dst &dest, double epsilon) const
        {
            ptic("HsMetric::ProjectSparse");
            
            size_t nRows = topLeftNumRows();

            if (!factorizedLaplacian.initialized)
            {
                // Assemble the cotan Laplacian
                std::vector<Triplet> triplets, triplets3x;
                H1::getTriplets(triplets, mesh, geom, epsilon);
                // Expand the matrix by 3x
                MatrixUtils::TripleTriplets(triplets, triplets3x);

                // Add constraint rows / cols for "simple" constraints included in Laplacian
                addSimpleConstraintTriplets(triplets3x);
                // Pre-factorize the cotan Laplacian
                Eigen::SparseMatrix<double> L(nRows, nRows);
                L.setFromTriplets(triplets3x.begin(), triplets3x.end());
                factorizedLaplacian.Compute(L);
            }

            // Multiply by L^{-1} once by solving Lx = b
            Eigen::VectorXd mid = factorizedLaplacian.Solve(gradientCol);

            if (!bvh)
            {
                throw std::runtime_error("Must have a BVH to use sparse approximation");
            }

            else
            {
                getBlockClusterTree()->MultiplyV3(mid, mid, BCTKernelType::FractionalOnly);
            }

            // Re-zero out Lagrange multipliers, since the first solve
            // will have left some junk in them
            for (size_t i = 3 * mesh->nVertices(); i < nRows; i++)
            {
                mid(i) = 0;
            }

            // Multiply by L^{-1} again by solving Lx = b
            dest = factorizedLaplacian.Solve(mid);
            
            ptoc("HsMetric::ProjectSparse");
        }

        template <typename Inverse, typename HsPtr>
        void GetSchurComplement(const HsPtr hs, SchurComplement &dest)
        {
            ptic("GetSchurComplement");
            
            size_t nVerts = hs->mesh->nVertices();
            size_t compNRows = 0;
            size_t bigNRows = hs->topLeftNumRows();

            // Figure out how many rows the constraint block is
            for (const ConstraintPack &c : hs->newtonConstraints)
            {
                compNRows += c.constraint->nRows();
            }
            if (compNRows == 0)
            {
                throw std::runtime_error("No constraints provided to Schur complement.");
            }

            dest.C.setZero(compNRows, bigNRows);
            size_t curRow = 0;

            // Fill in the constraint block by getting the entries for each constraint
            // while incrementing the rows
            for (const ConstraintPack &c : hs->newtonConstraints)
            {
                c.constraint->addEntries(dest.C, hs->mesh, hs->geom, curRow);
                curRow += c.constraint->nRows();
            }

            // https://en.wikipedia.org/wiki/Schur_complement
            // We want to compute (M/A) = D - C A^{-1} B.
            // In our case, D = 0, and B = C^T, so this is -C A^{-1} C^T.
            // Unfortunately this means we have to apply A^{-1} once to each column of C^T,
            // which could get expensive if we have too many constraints.

            // First allocate some space for a single column
            Eigen::VectorXd curCol;
            curCol.setZero(bigNRows);
            // And some space for A^{-1} C^T
            dest.Ainv_CT.setZero(bigNRows, compNRows);

            // For each column, copy it into curCol, and do the solve for A^{-1}
            for (size_t r = 0; r < compNRows; r++)
            {
                // Copy the row of C into the column
                for (size_t i = 0; i < 3 * nVerts; i++)
                {
                    curCol(i) = dest.C(r, i);
                }
                std::cout << "  Applying metric inverse to compute Schur complement row " << (r + 1) << "..." << std::endl;
                Inverse::Apply(*hs, curCol, curCol);
                // Copy the column into the column of A^{-1} C^T
                for (size_t i = 0; i < bigNRows; i++)
                {
                    dest.Ainv_CT(i, r) = curCol(i);
                }
            }

            // Now we've multiplied A^{-1} C^T, so just multiply this with C and negate it
            dest.M_A = -dest.C * dest.Ainv_CT;
            
            ptoc("GetSchurComplement");
        }

        class SparseInverse
        {
        public:
            template <typename V, typename Dest>
            static void Apply(const HsMetric &hs, const V &gradient, Dest &dest)
            {
                hs.InvertMetric(gradient, dest);
            }
        };
    } // namespace Hs

} // namespace biorsurfaces
