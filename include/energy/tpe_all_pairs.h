#pragma once

#include "rsurface_types.h"
#include "surface_energy.h"
#include "helpers.h"
#include "optimized_bct_types.h"
#include "optimized_bct.h"
#include "derivative_assembler.h"


namespace biorsurfaces
{

    // 0-th order multipole approximation of tangent-point energy
    class TPEnergyAllPairs : public SurfaceEnergy
    {
    public:
        TPEnergyAllPairs( MeshPtr mesh_, GeomPtr geom_, mreal alpha_, mreal beta_,  mreal weight_ = 1.)
        {
            bvh = 0;
            alpha = alpha_;
            beta = beta_;
            weight = weight_;
            mesh = mesh_;
            geom = geom_;
            
            mreal intpart;
            use_int = (std::modf( alpha, &intpart) == 0.0) && (std::modf( beta/2, &intpart) == 0.0);

            Update();
        }

        ~TPEnergyAllPairs()
        {
            if (bvh) delete bvh;
        }
        
        // Returns the current value of the energy.
        virtual double Value();
        // Returns the current differential of the energy, stored in the given
        // V x 3 matrix, where each row holds the differential (a 3-vector) with
        // respect to the corresponding vertex.
        virtual void Differential(Eigen::MatrixXd &output);

        // Update the energy to reflect the current state of the mesh. This could
        // involve building a new BVH for Barnes-Hut energies, for instance.
        virtual void Update();
        
        // Get the exponents of this energy; only applies to tangent-point energies.
        virtual Vector2 GetExponents();

        // Get a pointer to the current BVH for this energy.
        // Return 0 if the energy doesn't use a BVH.
        virtual OptimizedClusterTree *GetBVH();

        // Return the separation parameter for this energy.
        // Return 0 if this energy doesn't do hierarchical approximation.
        virtual double GetTheta();
        
        OptimizedBlockClusterTree * GetBCT();
        
        bool use_int = false;
        
    private:
        OptimizedClusterTree* bvh = nullptr;
        
        mreal alpha = 6.;
        mreal beta  = 12.;
        
        template<typename T1, typename T2>
        mreal Energy( T1 alpha, T2 betahalf);
        
        template<typename T1, typename T2>
        mreal DEnergy(T1 alpha, T2 betahalf);
        
    }; // TPEnergyMultipole0

} // namespace biorsurfaces
