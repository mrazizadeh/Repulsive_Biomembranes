#include "sobolev/constraints/total_area.h"
#include "surface_derivatives.h"
#include "helpers.h"
#include <cmath>
namespace biorsurfaces
{
    namespace Constraints
    {
        TotalAreaConstraint::TotalAreaConstraint(const MeshPtr &mesh, const GeomPtr &geom)
        {
            ResetFunction(mesh, geom);
        }

        void TotalAreaConstraint::ResetFunction(const MeshPtr &mesh, const GeomPtr &geom)
        {   
            initValue = totalArea(geom, mesh)*targetValue;            //std::cout<<"######Target Value is#########: "<<targetValue<<std::endl;
        }

        size_t TotalAreaConstraint::nRows()
        {
            return 1;
        }

        // Derivative of total volume is the mean curvature normal at each vertex.
        void TotalAreaConstraint::addTriplets(std::vector<Triplet> &triplets, const MeshPtr &mesh, const GeomPtr &geom, int baseRow)
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

        void TotalAreaConstraint::addEntries(Eigen::MatrixXd &M, const MeshPtr &mesh, const GeomPtr &geom, int baseRow)
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

        void TotalAreaConstraint::addErrorValues(Eigen::VectorXd &V, const MeshPtr &mesh, const GeomPtr &geom, int baseRow)
        {
            double current = totalArea(geom, mesh);
            V(baseRow) = (current - initValue);
        }

        double TotalAreaConstraint::getTargetValue()
        {
            return initValue;
        }

        void TotalAreaConstraint::setTargetValue(double targetValue_)
        {            
            //std::cout<<"######Target Value is changed from #########: "<<targetValue<<std::endl;

            targetValue = targetValue_;
            //std::cout<<"######TO: #########: "<<targetValue<<std::endl;

        }

        void TotalAreaConstraint::incrementTargetValue(double incr)
        {
            initValue += incr;
        }

    } // namespace Constraints
} // namespace biorsurfaces