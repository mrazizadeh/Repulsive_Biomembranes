#pragma once

#include "rsurface_types.h"
#include "matrix_utils.h"

namespace biorsurfaces
{
    namespace Constraints
    {
        class ConstraintBase
        {
        public:
            virtual ~ConstraintBase() {}
            virtual void ResetFunction(const MeshPtr &mesh, const GeomPtr &geom) = 0;
            virtual void addTriplets(std::vector<Triplet> &triplets, const MeshPtr &mesh, const GeomPtr &geom, int baseRow) = 0;
            virtual void addEntries(Eigen::MatrixXd &M, const MeshPtr &mesh, const GeomPtr &geom, int baseRow) = 0;
            virtual void addErrorValues(Eigen::VectorXd &V, const MeshPtr &mesh, const GeomPtr &geom, int baseRow) = 0;
            virtual size_t nRows() = 0;
        };

        void addEntriesToSymmetric(ConstraintBase &cs, Eigen::MatrixXd &M, const MeshPtr &mesh, const GeomPtr &geom, int baseRow);
        void addTripletsToSymmetric(ConstraintBase &cs, std::vector<Triplet> &triplets, const MeshPtr &mesh, const GeomPtr &geom, int baseRow);

        class SaddleMatrixConstraint : public ConstraintBase
        {
        public:
            virtual double getTargetValue() = 0;
            virtual void setTargetValue(double targetValue_)=0;
            virtual void incrementTargetValue(double incr) = 0;
        };

        class SimpleProjectorConstraint : public ConstraintBase
        {
        public:
            virtual void ProjectConstraint(MeshPtr &mesh, GeomPtr &geom) = 0;
        };

    } // namespace Constraints

    struct ConstraintPack
    {
        Constraints::SaddleMatrixConstraint *constraint;
        double stepSize;
        long iterationsLeft;
    };

} // namespace biorsurfaces