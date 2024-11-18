#pragma once

#include "sobolev/constraints.h"

namespace biorsurfaces
{
    namespace Constraints
    {
        class ScalingConstraint : public SaddleMatrixConstraint
        {
        public:
            virtual void ResetFunction(const MeshPtr &mesh, const GeomPtr &geom);
            virtual void addTriplets(std::vector<Triplet> &triplets, const MeshPtr &mesh, const GeomPtr &geom, int baseRow);
            virtual void addEntries(Eigen::MatrixXd &M, const MeshPtr &mesh, const GeomPtr &geom, int baseRow);
            virtual void addValue(Eigen::VectorXd &V, const MeshPtr &mesh, const GeomPtr &geom, int baseRow);
            virtual double getTargetValue();
            virtual void setTargetValue(double targetValue_);
            virtual void incrementTargetValue(double incr);
            virtual size_t nRows();
        };
    } // namespace Constraints
} // namespace biorsurfaces