#include "sobolev/constraints/reduced_volume.h"
#include "surface_derivatives.h"
#include "helpers.h"
#include <cmath>
namespace biorsurfaces
{
    namespace Constraints
    {
        ReducedVolume::ReducedVolume(const MeshPtr &mesh, const GeomPtr &geom)
        {
            ResetFunction(mesh, geom);
        }

        void ReducedVolume::ResetFunction(const MeshPtr &mesh, const GeomPtr &geom)
        {
            volume = totalVolume(geom, mesh);
            area = totalArea(geom, mesh);
                // std::cout<<"######Target Value is#########: "<<targetValue<<std::endl;
                // double targetValue = 0.642;
                initValue = pow(volume / targetValue, 2.0 / 3.0) * 4.83597586205;
            // initValue = targetValue*totalArea(geom, mesh);
            // std::cout<<"######Init Values is#########: "<<initValue<<std::endl;
        }

        size_t ReducedVolume::nRows()
        {
            return 1;
        }

        // Derivative of total volume is the mean curvature normal at each vertex.
        void ReducedVolume::addTriplets(std::vector<Triplet> &triplets, const MeshPtr &mesh, const GeomPtr &geom, int baseRow)
        {
            VertexIndices indices = mesh->getVertexIndices();
            for (GCVertex v : mesh->vertices())
            {
                Vector3 normal_v = SurfaceDerivs::meanCurvatureNormal(v, geom);
                size_t i3 = 3 * indices[v];
                // Fill constraint row with mean curvature normal in each vertex's 3 entries
                triplets.push_back(Triplet(baseRow, i3, normal_v.x));
                triplets.push_back(Triplet(baseRow, i3 + 1, normal_v.y));
                triplets.push_back(Triplet(baseRow, i3 + 2, normal_v.z));
            }
        }

        void ReducedVolume::addEntries(Eigen::MatrixXd &M, const MeshPtr &mesh, const GeomPtr &geom, int baseRow)
        {
            VertexIndices indices = mesh->getVertexIndices();
            for (GCVertex v : mesh->vertices())
            {
                Vector3 normal_v = SurfaceDerivs::meanCurvatureNormal(v, geom);
                size_t i3 = 3 * indices[v];
                // Fill constraint row with mean curvature normal in each vertex's 3 entries
                M(baseRow, i3) = normal_v.x;
                M(baseRow, i3 + 1) = normal_v.y;
                M(baseRow, i3 + 2) = normal_v.z;
            }
        }

        void ReducedVolume::addErrorValues(Eigen::VectorXd &V, const MeshPtr &mesh, const GeomPtr &geom, int baseRow)
        {
            double current = totalArea(geom, mesh);
            V(baseRow) = (current - initValue);
        }

        double ReducedVolume::getTargetValue()
        {
            return volume / ((4 / 3.0) * geometrycentral::PI * std::sqrt(area / 4 / geometrycentral::PI) * std::sqrt(area / 4 / geometrycentral::PI) * std::sqrt(area / 4 / geometrycentral::PI));
        }

        void ReducedVolume::setTargetValue(double targetValue_)
        {
            // std::cout<<"######Target Value is changed from #########: "<<targetValue<<std::endl;

            targetValue = targetValue_;
            // std::cout<<"######TO: #########: "<<targetValue<<std::endl;
        }

        void ReducedVolume::incrementTargetValue(double incr)
        {
            initValue += incr;
        }

    } // namespace Constraints
} // namespace biorsurfaces