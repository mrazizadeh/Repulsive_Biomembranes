#pragma once

#include "rsurface_types.h"
#include "optimized_bct.h"
#include "metric_term.h"
#include "sobolev/hs.h"

class BCTMatrixReplacement;
using Eigen::SparseMatrix;

// Adapted from https://eigen.tuxfamily.org/dox/group__MatrixfreeSolverExample.html

namespace Eigen
{
    namespace internal
    {
        // MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
        template <>
        struct traits<BCTMatrixReplacement> : public Eigen::internal::traits<Eigen::SparseMatrix<double>>
        {
        };
    } // namespace internal
} // namespace Eigen

class BCTMatrixReplacement : public Eigen::EigenBase<BCTMatrixReplacement>
{
public:
    // Required typedefs, constants, and method:
    typedef double Scalar;
    typedef double RealScalar;
    typedef int StorageIndex;
    enum
    {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
    };

    Eigen::Index rows() const
    {
        return 3 * hs->mesh->nVertices() + C->rows();
    }

    Eigen::Index outerSize() const
    {
        return 3 * hs->mesh->nVertices() + C->rows();
    }

    Eigen::Index innerSize() const
    {
        return 3 * hs->mesh->nVertices() + C->rows();
    }

    Eigen::Index cols() const
    {
        return 3 * hs->mesh->nVertices() + C->rows();
    }

    template <typename Rhs>
    Eigen::Product<BCTMatrixReplacement, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs> &x) const
    {
        return Eigen::Product<BCTMatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
    }

    // Custom API:
    BCTMatrixReplacement() {}

    void addConstraintBlock(const Eigen::SparseMatrix<double> &C_)
    {
        C = &C_;
    }

    void addMetric(const biorsurfaces::Hs::HsMetric *hs_)
    {
        hs = hs_;
    }

    const biorsurfaces::Hs::HsMetric *getHs() const
    {
        return hs;
    }

    const Eigen::SparseMatrix<double> &getConstraintBlock() const
    {
        return *C;
    }

private:
    const Eigen::SparseMatrix<double> *C;
    const biorsurfaces::Hs::HsMetric *hs;
};

namespace biorsurfaces
{

    namespace Hs
    {
        class SparseHsPreconditioner
        {
            typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

        public:
            typedef typename Vector::StorageIndex StorageIndex;
            enum
            {
                ColsAtCompileTime = Eigen::Dynamic,
                MaxColsAtCompileTime = Eigen::Dynamic
            };
            SparseHsPreconditioner() {}

            template <typename MatrixType>
            explicit SparseHsPreconditioner(const MatrixType &fracL)
            {
                compute(fracL);
            }

            Eigen::Index rows() const { return hs->topLeftNumRows(); }
            Eigen::Index cols() const { return hs->topLeftNumRows(); }

            template <typename MatrixType>
            SparseHsPreconditioner &analyzePattern(const MatrixType &) { return *this; }

            template <typename MatrixType>
            SparseHsPreconditioner &factorize(const MatrixType &) { return *this; }

            template <typename MatrixType>
            SparseHsPreconditioner &compute(const MatrixType &fracL)
            {
                hs = fracL.getHs();
                return *this;
            }

            /** \internal */
            template <typename Rhs, typename Dest>
            void _solve_impl(const Rhs &b, Dest &x) const
            {
                std::cout << "  * GMRES iteration " << (count++) << "...\r" << std::flush;
                x = hs->InvertSparseForIterative(b);
            }

            template <typename Rhs>
            inline const Eigen::Solve<SparseHsPreconditioner, Rhs>
            solve(const Eigen::MatrixBase<Rhs> &b) const
            {
                return Eigen::Solve<SparseHsPreconditioner, Rhs>(*this, b.derived());
            }

            const Hs::HsMetric *hs;
            const std::vector<ConstraintPack> schurConstraints;
            SchurComplement schur;

            mutable size_t count = 0;

            Eigen::ComputationInfo info() { return Eigen::Success; }
        };
    } // namespace Hs
} // namespace biorsurfaces

namespace Eigen
{
    namespace internal
    {

        template <typename Rhs>
        struct generic_product_impl<BCTMatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
            : generic_product_impl_base<BCTMatrixReplacement, Rhs, generic_product_impl<BCTMatrixReplacement, Rhs>>
        {
            typedef typename Product<BCTMatrixReplacement, Rhs>::Scalar Scalar;

            template <typename Dest>
            static void scaleAndAddTo(Dest &dst, const BCTMatrixReplacement &lhs, const Rhs &rhs, const Scalar &alpha)
            {
                // This method should implement "dst += alpha * lhs * rhs" inplace,
                // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
                assert(alpha == Scalar(1) && "scaling is not implemented");
                EIGEN_ONLY_USED_FOR_DEBUG(alpha);

                const biorsurfaces::Hs::HsMetric* hs = lhs.getHs();

                Eigen::VectorXd product(3 * hs->mesh->nVertices() + lhs.getConstraintBlock().rows());
                product.setZero();

                Eigen::VectorXd rhsCopy = rhs;

                for (biorsurfaces::MetricTerm* term : hs->getMetricTerms())
                {
                    term->MultiplyAdd(rhsCopy, product);
                }
                biorsurfaces::MultiplyConstraintBlock(hs->mesh, rhsCopy, product, lhs.getConstraintBlock(), true);

                dst += product;
            }
        };

    } // namespace internal
} // namespace Eigen
