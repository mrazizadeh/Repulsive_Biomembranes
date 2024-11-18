#pragma once

#include "rsurface_types.h"
#include "energy/tpe_kernel.h"
#include "spatial/bvh_6d.h"

namespace biorsurfaces
{
    class CoulombEnergy : public SurfaceEnergy
    {
    public:
        CoulombEnergy(TPEKernel *kernel_, double theta_);
        ~CoulombEnergy();
        virtual double Value();
        virtual void Differential(Eigen::MatrixXd &output);
        virtual void Update();
        virtual MeshPtr GetMesh();
        virtual GeomPtr GetGeom();
        virtual Vector2 GetExponents();
        virtual OptimizedClusterTree *GetBVH();
        virtual double GetTheta();

    private:
        TPEKernel *kernel;
        OptimizedClusterTree *root;
        double theta;
        double energyAtVertex(OptimizedClusterTree *node, GCVertex v);
        Vector3 gradientAtVertex(OptimizedClusterTree *node, GCVertex v);
    };
} // namespace biorsurfaces