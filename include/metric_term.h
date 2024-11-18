#pragma once

#include "rsurface_types.h"
#include "optimized_bct.h"

namespace biorsurfaces
{
    class MetricTerm
    {
        public:
        virtual ~MetricTerm() {}
        // Multiply the metric term with vec, and add the product to result.
        virtual void MultiplyAdd(Eigen::VectorXd &vec, Eigen::VectorXd &result) const = 0;
    };

    class BCTMetricTerm : public MetricTerm
    {
        public:
        BCTMetricTerm(std::shared_ptr<OptimizedBlockClusterTree> bct_)
        {
            bct = bct_;
        }
        
        virtual void MultiplyAdd(Eigen::VectorXd &vec, Eigen::VectorXd &result) const
        {
            bct->MultiplyV3(vec, result, BCTKernelType::HighOrder, true);
            bct->MultiplyV3(vec, result, BCTKernelType::LowOrder, true);
        }

        private:
        std::shared_ptr<OptimizedBlockClusterTree> bct;
    };

    class BiLaplacianMetricTerm : public MetricTerm
    {
        public:
        BiLaplacianMetricTerm(MeshPtr &mesh, GeomPtr &geom);
        virtual void MultiplyAdd(Eigen::VectorXd &vec, Eigen::VectorXd &result) const;

        private:
        size_t nMultiplyRows;
        Eigen::SparseMatrix<double> biLaplacian;
    };
}

