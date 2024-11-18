#pragma once

#include "energy/tpe_kernel.h"
#include "spatial/bvh_6d.h"

namespace biorsurfaces
{

    // Evaluates energy and differential using Barnes-Hut with a BVH.
    class BarnesHutTPEnergy6D : public SurfaceEnergy
    {
    public:
        BarnesHutTPEnergy6D(TPEKernel *kernel_, double theta_);
        ~BarnesHutTPEnergy6D();
        virtual double Value();
        virtual void Differential(Eigen::MatrixXd &output);
        virtual void Update();
        virtual MeshPtr GetMesh();
        virtual GeomPtr GetGeom();
        virtual Vector2 GetExponents();
        virtual OptimizedClusterTree *GetBVH();
        virtual double GetTheta();
        geometrycentral::surface::FaceData<double> energyPerFace;

        bool disableNearField;

    private:
        TPEKernel *kernel;
        BVHNode6D *root;

        double theta;
        double computeEnergyOfFace(GCFace face, BVHNode6D *bvhRoot);
        void accumulateTPEGradient(Eigen::MatrixXd &gradients, BVHNode6D *node, GCFace face1,
                                   surface::VertexData<size_t> &indices);
    };

} // namespace biorsurfaces