#pragma once

#include <iostream>
#include <vector>
#include <sstream>
#include "rsurface_types.h"
#include "implicit/simple_surfaces.h"

namespace biorsurfaces
{
    bool endsWith(std::string const &fullString, std::string const &ending);

    namespace scene
    {
        enum class ConstraintType
        {
            Barycenter,
            TotalArea,
            TotalVolume,
            BoundaryPins,
            VertexPins,
            BoundaryNormals,
            VertexNormals,
            ReducedVolume,
            ReducedMC 
        };

        enum class PotentialType
        {
            SquaredError,
            BoundaryLength,
            BoundaryCurvature,
            Area,
            Volume,
            SoftAreaConstraint,
            SoftVolumeConstraint,
            Helfrich,
            HelfrichADE,
            SoftTotalMeanCurvatureConstraint,
            SoftReducedVolumeConstraint
        };

        struct PotentialData
        {
            PotentialType type;
            double weight;
            double targetValue;
        };

        struct ConstraintData
        {
            ConstraintType type;
            double targetMultiplier;
            long numIterations;
            double targetAddition;
        };

        struct ObstacleData
        {
            std::string obstacleName;
            double weight;
            bool recenter = false;
            bool asPointCloud = false;
        };

        enum class ImplicitType
        {
            Sphere,
            Cylinder,
            Torus,
            Plane
        };

        struct ImplicitBarrierData
        {
            ImplicitType type;
            std::vector<double> parameters;
            bool repel = true;
            double power = 2;
            double weight = 1;
        };

        struct SceneData
        {
            std::string meshName;
            double alpha;
            double beta;
            bool allowBarycenterShift = false;
            std::vector<ObstacleData> obstacles;
            std::vector<ConstraintData> constraints;
            std::vector<PotentialData> potentials;
            std::vector<VertexPinData> vertexPins;
            std::vector<size_t> vertexNormals;
            std::vector<ImplicitBarrierData> implicitBarriers;
            int iterationLimit = 5000;
            long realTimeLimit = 0;
            std::string performanceLogFile = "performance.csv";
            GradientMethod defaultMethod = GradientMethod::HsProjectedIterative;
            bool disableNearField = false;
            bool autoComputeVolumeTarget = false;
            double autoVolumeTargetRatio = 1;
        };

        template <class Container>
        void splitString(const std::string &str, Container &cont, char delim = ' ')
        {
            std::stringstream ss(str);
            std::string token;
            while (std::getline(ss, token, delim))
            {
                cont.push_back(token);
            }
        }
        SceneData parseScene(std::string filename);

    } // namespace scene

    

} // namespace biorsurfaces
