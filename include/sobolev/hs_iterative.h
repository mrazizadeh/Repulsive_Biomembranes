#pragma once

#include "sobolev/hs.h"
#include "bct_matrix_replacement.h"
#include "metric_term.h"

#include <unsupported/Eigen/IterativeSolvers>

namespace biorsurfaces
{
    namespace Hs
    {
        void ProjectConstrainedHsIterativeMat(Hs::HsMetric &hs, Eigen::MatrixXd &gradient, Eigen::MatrixXd &dest);
        void ProjectConstrainedHsIterative(Hs::HsMetric &hs, Eigen::VectorXd &gradient, Eigen::VectorXd &dest);

        template <typename V, typename Dest>
        void ProjectUnconstrainedHsIterative(const Hs::HsMetric &hs, const V &gradient, Dest &dest,
                                             Eigen::SparseMatrix<double> &constraintBlock)
        {
            ptic("ProjectUnconstrainedHsIterative");
            
            BCTPtr bct = hs.getBlockClusterTree();

            BCTMatrixReplacement fracL;
            fracL.addConstraintBlock(constraintBlock);
            fracL.addMetric(&hs);

            // Eigen::ConjugateGradient<BCTMatrixReplacement, Eigen::Lower | Eigen::Upper, SparseHsPreconditioner> cg;
            Eigen::GMRES<BCTMatrixReplacement, SparseHsPreconditioner> cg;
            cg.compute(fracL);
            
            Eigen::VectorXd temp;
            temp.setZero(gradient.rows());
            cg.setTolerance(1e-4);
            temp = cg.solveWithGuess(gradient, temp);
            std::cout << "  * GMRES converged in " << cg.iterations() << " iterations, final residual = " << cg.error() << std::endl;

            dest = temp;
            
            ptoc("ProjectUnconstrainedHsIterative");
        }

        template <typename V, typename Dest>
        void ProjectUnconstrainedHsIterative(const Hs::HsMetric &hs, const V &gradient, Dest &dest)
        {
            Eigen::SparseMatrix<double> constraintBlock = hs.GetConstraintBlock(false);
            ProjectUnconstrainedHsIterative(hs, gradient, dest, constraintBlock);
        }

        class IterativeInverse
        {
        public:
            template <typename V, typename Dest>
            static void Apply(const HsMetric &hs, const V &gradient, Dest &dest)
            {
                ProjectUnconstrainedHsIterative(hs, gradient, dest);
            }
        };
    } // namespace Hs
} // namespace biorsurfaces
