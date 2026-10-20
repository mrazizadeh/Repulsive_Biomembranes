#pragma once

#include "hs.h"
#include "hard_exclusion/bridge.h"
#include "hard_exclusion/constraint_projection.h"

namespace biorsurfaces
{

    namespace Hs
    {
        template <typename Inverse>
        void ProjectViaSchurV(const HsMetric &hs, Eigen::VectorXd &curCol, Eigen::VectorXd &dest);

        template <typename Inverse>
        void ProjectViaSchur(const HsMetric &hs, Eigen::MatrixXd &gradient, Eigen::MatrixXd &dest);

        template <typename Inverse>
        void ProjectSchurConstraints(const HsMetric &hs, int newtonSteps);

        template <typename Inverse>
        void ProjectViaSchurV(const HsMetric &hs, Eigen::VectorXd &curCol, Eigen::VectorXd &dest)
        {
            size_t nVerts = hs.mesh->nVertices();
            // Invert the "saddle matrix" now:
            // the block of M^{-1} we want is A^{-1} + A^{-1} C^T (M/A)^{-1} C A^{-1}
            Eigen::VectorXd Ainv_g = curCol;
            std::cout << "  Applying metric inverse for gradient..." << std::endl;
            Inverse::Apply(hs, Ainv_g, Ainv_g);
            std::cout << "  Applied." << std::endl;

            // Now we compute the correction
            // Start from hsGradient = A^{-1} x
            Eigen::VectorXd C_Ai_x = hs.Schur<Inverse>().C * Ainv_g;
            Eigen::VectorXd MAi_C_Ai_x;
            MAi_C_Ai_x.setZero(C_Ai_x.rows());
            MatrixUtils::SolveDenseSystem(hs.Schur<Inverse>().M_A, C_Ai_x, MAi_C_Ai_x);
            // Use the cached A^{-1} C^T to do this last step
            std::cout << "  Applying metric inverse for Schur complement orthogonalization..." << std::endl;
            dest = Ainv_g + hs.Schur<Inverse>().Ainv_CT * MAi_C_Ai_x;
        }

        template <typename Inverse>
        void ProjectViaSchur(const HsMetric &hs, Eigen::MatrixXd &gradient, Eigen::MatrixXd &dest)
        {
            size_t nVerts = hs.mesh->nVertices();

            if (hs.newtonConstraints.size() > 0)
            {
                Eigen::VectorXd temp;
                temp.setZero(hs.Schur<Inverse>().C.cols());
                MatrixUtils::MatrixIntoColumn(gradient, temp);

                ProjectViaSchurV<Inverse>(hs, temp, temp);

                MatrixUtils::ColumnIntoMatrix(temp, dest);
            }
            else
            {
                hs.InvertMetricMat(gradient, dest);
            }
        }

        template <typename Inverse>
        void ProjectSchurConstraints(const HsMetric &hs, int newtonSteps)
        {
            ptic("ProjectSchurConstraints");
            
            size_t nRows = hs.Schur<Inverse>().M_A.rows();
            Eigen::VectorXd vals(nRows);

            const auto correctOnce = [&](int nIters)
            {
                vals.setZero();
                int curRow = 0;
                // Fill right-hand side with error values
                for (const ConstraintPack &c : hs.newtonConstraints)
                {
                    c.constraint->addErrorValues(vals, hs.mesh, hs.geom, curRow);
                    curRow += c.constraint->nRows();
                }

                double constraintError = vals.lpNorm<Eigen::Infinity>();
                std::cout << "  * Constraint error after " << nIters << " iterations = " << constraintError << std::endl;
                if (!HardExclusion::Enabled() && nIters > 0 && constraintError < 1e-2)
                {
                    return false;
                }

                // In this case we want the block of the inverse that multiplies the bottom block
                // -A^{-1} B (M/A)^{-1}, where B = C^T
                // Apply (M/A) inverse first
                MatrixUtils::SolveDenseSystem(hs.Schur<Inverse>().M_A, vals, vals);
                // Apply cached A{^-1} C^T
                Eigen::VectorXd correction = hs.Schur<Inverse>().Ainv_CT * vals;
                if(HardExclusion::Enabled() && !correction.allFinite())
                    throw std::runtime_error("nonfinite_constraint_correction");

                // Independently guard Newton correction; guarding descent alone is insufficient.
                if(HardExclusion::Enabled()) {
                    auto before=HardExclusion::Capture(hs.mesh,hs.geom);
                    hard_exclusion::Positions proposed=before;
                    for(int i=0;i<proposed.rows();++i) for(int j=0;j<3;++j) proposed(i,j)+=correction(3*i+j);
                    correction*=HardExclusion::Limit(before,proposed,"constraint_correction");
                }
                // Apply the correction to the vertex positions
                VertexIndices verts = hs.mesh->getVertexIndices();
                for (GCVertex v : hs.mesh->vertices())
                {
                    int base = 3 * verts[v];
                    Vector3 vertCorr{correction(base), correction(base + 1), correction(base + 2)};
                    hs.geom->inputVertexPositions[v] += vertCorr;
                }
                if(HardExclusion::Enabled()) hs.geom->refreshQuantities();

                vals.setZero();
                curRow = 0;
                // Fill right-hand side with error values
                for (const ConstraintPack &c : hs.newtonConstraints)
                {
                    c.constraint->addErrorValues(vals, hs.mesh, hs.geom, curRow);
                    curRow += c.constraint->nRows();
                }
                double corrError = vals.lpNorm<Eigen::Infinity>();
                std::cout << "  * Corrected error " << constraintError << " -> " << corrError << std::endl;
                return true;
            };
            if(HardExclusion::Enabled()) {
                // The legacy caller requests only one correction. Repeat the
                // same cached Schur correction for guarded runs until the
                // declared geometric A/V tolerance holds, or stop explicitly.
                // No relaxation of the full-step tolerance or energy gate.
                const int corrections=HardExclusion::CorrectConstraintsToTolerance(
                    [&]() { return HardExclusion::ConstraintRelativeError(hs.mesh,hs.geom); },
                    [&](int i) { correctOnce(i); },
                    HardExclusion::ConstraintRelativeTolerance());
                std::cout << "  * Relative A/V error after " << corrections << " guarded corrections = "
                          << HardExclusion::ConstraintRelativeError(hs.mesh,hs.geom) << std::endl;
            } else {
                // Preserve the original count and absolute stopping criterion
                // exactly when the independent guard is disabled.
                for(int i=0;i<newtonSteps;++i) if(!correctOnce(i)) break;
            }
            
            ptoc("ProjectSchurConstraints");
        }
    } // namespace Hs

} // namespace biorsurfaces
