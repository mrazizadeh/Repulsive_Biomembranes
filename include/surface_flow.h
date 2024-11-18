#pragma once

#include <memory>

#include "rsurface_types.h"
#include "energy/tp_obstacle_barnes_hut_0.h"
#include "sobolev/constraints.h"
#include "sobolev/all_constraints.h"
#include "surface_energy.h"
#include "line_search.h"
#include "sobolev/hs_ncg.h"
#include "sobolev/lbfgs.h"
#include "profiler.h"

namespace biorsurfaces
{
    class SurfaceFlow
    {
    public:
        SurfaceFlow(SurfaceEnergy *energy_);
        void AddAdditionalEnergy(SurfaceEnergy *extraEnergy);
        void AddObstacleEnergy(SurfaceEnergy *obsEnergy);

        void StepL2Unconstrained();
        void StepL2Projected();
        void StepProjectedGradientExact();
        void StepProjectedGradient();
        void StepProjectedGradientIterative();
        void StepH1LBFGS();
        void StepBQN();
        void StepH1ProjGrad();
        void StepAQP(double invKappa);
        void StepH2Projected();

        SurfaceEnergy *BaseEnergy();

        void RecenterMesh();
        void ResetAllConstraints();
        //INACTIVATED BY MEGHDAD-Jan 3, 2023
        void ResetAllPotentials();

        void AssembleGradients(Eigen::MatrixXd &dest);
        std::unique_ptr<Hs::HsMetric> GetHsMetric();

        void UpdateEnergies();
        double evaluateEnergy();

        template <typename Constraint>
        Constraint *addSchurConstraint(MeshPtr &mesh, GeomPtr &geom, double multiplier, long iterations, double add = 0)
        {
            Constraint *c = new Constraint(mesh, geom);
            double stepSize = 0;
            if (iterations > 0)
            {
                double initVal = c->getTargetValue();
                double target = multiplier * initVal + add;
                double change = target - initVal;
                std::cout << "Target value for constraint is " << target << " (initial = " << initVal << ")" << std::endl;
                stepSize = change / iterations;
            }
            schurConstraints.push_back(ConstraintPack{c, stepSize, iterations});
            return c;
        }

        template <typename Constraint>
        Constraint *addSimpleConstraint(MeshPtr &mesh, GeomPtr &geom)
        {
            Constraint *c = new Constraint(mesh, geom);
            simpleConstraints.push_back(c);
            return c;
        }

        template <typename Constraint>
        inline ConstraintPack& findPackOfSchurConstraintType()
        {
            for (ConstraintPack &c : schurConstraints)
            {
                if (dynamic_cast<Constraint*>(c.constraint))
                {
                    return c;
                }
            }
            throw std::runtime_error("No Schur constraint of the given type found.");
        }

        inline void retargetSchurConstraint(ConstraintPack &c, double newValue)
        {
            double curTarget = c.constraint->getTargetValue();
            double diff = newValue - curTarget;
            c.stepSize = diff / c.iterationsLeft;
        }

        template <typename Constraint>
        inline void retargetSchurConstraintOfType(double newValue)
        {
            retargetSchurConstraint(findPackOfSchurConstraintType<Constraint>(), newValue);
        }

        bool allowBarycenterShift;
        bool verticesMutated;
        bool disableNearField;

        // if this value is positive, the flow will not
        // take steps larger than the given value
        double maxStepSize = -1.;
        std::vector<SurfaceEnergy *> energies;
        std::vector<ConstraintPack> schurConstraints;


    private:
        //std::vector<SurfaceEnergy *> energies;
        MeshPtr mesh;
        GeomPtr geom;
        Eigen::MatrixXd prevPositions1;
        Eigen::MatrixXd prevPositions2;
        unsigned int stepCount;
        Vector3 origBarycenter;
        Constraints::BarycenterComponentsConstraint *secretBarycenter;
        LBFGSOptimizer* lbfgs;
        SurfaceEnergy* obstacleEnergy;
        std::vector<Constraints::SimpleProjectorConstraint *> simpleConstraints;

        size_t addConstraintTriplets(std::vector<Triplet> &triplets, bool includeSchur);
        
        void prefactorConstrainedLaplacian(SparseFactorization &factored, bool includeSchur);
        void prefactorConstrainedLaplacian(Eigen::SparseMatrix<double> &L, SparseFactorization &factored, bool includeSchur);

        inline void incrementSchurConstraints()
        {
            for (ConstraintPack &c : schurConstraints)
            {
                if (c.iterationsLeft > 0)
                {
                    c.iterationsLeft--;
                    c.constraint->incrementTargetValue(c.stepSize);
                }
            }
        }

        double bqn_B;

        double prevStep;
    };
} // namespace biorsurfaces
