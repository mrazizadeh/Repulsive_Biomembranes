#pragma once

#include "rsurface_types.h"
#include "surface_energy.h"
#include "helpers.h"
#include "bct_constructors.h"
#include "optimized_bct_types.h"
#include "optimized_bct.h"
#include "derivative_assembler.h"


namespace biorsurfaces
{

    // 0-th order multipole approximation of tangent-point energy
    class TPObstacleAllPairs : public SurfaceEnergy
    {
    public:
        TPObstacleAllPairs( MeshPtr mesh_, GeomPtr geom_, SurfaceEnergy *bvhSharedFrom_, MeshPtr &obsMesh, GeomPtr &obsGeom, mreal alpha_, mreal beta_, mreal weight_  = 1.)
        {
            mesh = mesh_;
            geom = geom_;
            bvh = 0;
            bvhSharedFrom = bvhSharedFrom_;
            o_bvh = CreateOptimizedBVH(obsMesh, obsGeom);
            alpha = alpha_;
            beta = beta_;
            
            weight = weight_;
            
            mreal intpart;
            use_int = (std::modf( alpha, &intpart) == 0.0) && (std::modf( beta/2, &intpart) == 0.0);
            
            Update();
        }

        ~TPObstacleAllPairs()
        {
            if (o_bvh)
            {
                delete o_bvh;
            }
        }
        
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
        
        virtual void Update();
        OptimizedBlockClusterTree * GetBCT();
        
        bool use_int = false;
        
    private:
        SurfaceEnergy* bvhSharedFrom;
        OptimizedClusterTree* bvh = nullptr;
        OptimizedClusterTree* o_bvh = nullptr;
        
        mreal alpha = 6.;
        mreal beta  = 12.;
        
        template<typename T1, typename T2>
        mreal Energy( T1 alpha, T2 betahalf);
        
        template<typename T1, typename T2>
        mreal DEnergy(T1 alpha, T2 betahalf);
        
    }; // TPEnergyMultipole0

} // namespace biorsurfaces
