#pragma once

#include "rsurface_types.h"
#include "lbfgs.h"
#include "sobolev/sparse_factorization.h"
#include "sobolev/constraints.h"

namespace biorsurfaces
{
    class H1_LBFGS : public LBFGSOptimizer
    {
        public:
        H1_LBFGS(size_t memSize_, std::vector<Constraints::SimpleProjectorConstraint *> simpleConstraints_);

        virtual void ApplyInnerProduct(Eigen::VectorXd &input, Eigen::VectorXd &output);
        virtual void ApplyInverseInnerProduct(Eigen::VectorXd &input, Eigen::VectorXd &output);
        virtual void SetUpInnerProduct(MeshPtr &mesh, GeomPtr &geom);

        protected:
        Eigen::SparseMatrix<double> L;
        Eigen::VectorXd tempVector;
        std::vector<Constraints::SimpleProjectorConstraint *> simpleConstraints;
    };
}