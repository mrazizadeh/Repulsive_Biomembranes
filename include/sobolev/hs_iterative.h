#pragma once

#include "sobolev/hs.h"
#include "bct_matrix_replacement.h"
#include "metric_term.h"
#include "hard_exclusion/bridge.h"
#include <limits>

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
            if(HardExclusion::Enabled() && !gradient.allFinite())
                HardExclusion::RecordLinearSolve(int(Eigen::NumericalIssue),0,
                    std::numeric_limits<double>::quiet_NaN(),1e-4,false,"rhs");
            
            BCTPtr bct = hs.getBlockClusterTree();

            BCTMatrixReplacement fracL;
            fracL.addConstraintBlock(constraintBlock);
            fracL.addMetric(&hs);

            // Eigen::ConjugateGradient<BCTMatrixReplacement, Eigen::Lower | Eigen::Upper, SparseHsPreconditioner> cg;
            Eigen::GMRES<BCTMatrixReplacement, SparseHsPreconditioner> cg;
            cg.setTolerance(1e-4); // preserve the inherited solver tolerance
            cg.compute(fracL);
            if(HardExclusion::Enabled() && cg.info()!=Eigen::Success)
                HardExclusion::RecordLinearSolve(int(cg.info()),0,
                    std::numeric_limits<double>::quiet_NaN(),cg.tolerance(),false,"setup");
            
            Eigen::VectorXd temp;
            temp.setZero(gradient.rows());
            temp = cg.solveWithGuess(gradient, temp);
            std::cout << "  * GMRES finished: Eigen status = " << int(cg.info()) << ", iterations = " << cg.iterations()
                      << ", reported residual = " << cg.error() << ", tolerance = " << cg.tolerance() << std::endl;
            HardExclusion::RecordLinearSolve(int(cg.info()),cg.iterations(),cg.error(),cg.tolerance(),temp.allFinite(),"solve");

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
