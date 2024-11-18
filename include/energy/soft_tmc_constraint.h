#pragma once

#include "rsurface_types.h"
#include "surface_energy.h"

namespace biorsurfaces
{
    class SoftTotalMeanCurvatureConstraint : public SurfaceEnergy
    {
    public:
        SoftTotalMeanCurvatureConstraint(MeshPtr mesh_, GeomPtr geom_, double weight_, double targetValue_);

        virtual ~SoftTotalMeanCurvatureConstraint() {}

        // Returns the current value of the energy.
        virtual double Value();

        // Returns the current differential of the energy, stored in the given
        // V x 3 matrix, where each row holds the differential (a 3-vector) with
        // respect to the corresponding vertex.
        virtual void Differential(Eigen::MatrixXd &output);

        // Get the exponents of this energy; only applies to tangent-point energies.
        virtual Vector2 GetExponents();

        // Get a pointer to the current BVH for this energy.
        // Return 0 if the energy doesn't use a BVH.
        virtual OptimizedClusterTree *GetBVH();

        // Return the separation parameter for this energy.
        // Return 0 if this energy doesn't do hierarchical approximation.
        virtual double GetTheta();
        virtual void ResetTargetsTMC(double newTargetValue,double newWeight);
    private:
        double initialArea;
        double targetValue;
    };
}