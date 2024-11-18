// Note: The majority of these functions are derived from the mem3dg package.Mem3DG is developed by Cuncheng Zhu, Christopher T. Lee, with contributions from others. 
// Development of Mem3DG is funded in part by AFOSR MURI FA9550-18-1-0051, and a Hartwell Foundation Postdoctoral Fellowship.
// The original source code can be found at https://github.com/RangamaniLabUCSD/Mem3DG
// The paper can be found at https://www.sciencedirect.com/science/article/pii/S2667074722000192

#include "energy/helfrich_energy_ADE.h"
#include "matrix_utils.h"
#include <Eigen/Core>

namespace biorsurfaces
{
     /**
     * @brief Computes the vector from a halfedge to its next vertex.
     * 
     * @param he The halfedge from which to compute the vector.
     * @param vpg The vertex position geometry containing vertex positions.
     * @return Vector3 The vector from the current vertex to the next vertex.
     */
    inline Vector3 vecFromHalfedge(const GCHalfedge he, const geometrycentral::surface::VertexPositionGeometry &vpg)
    {
        return vpg.inputVertexPositions[he.next().vertex()] -
               vpg.inputVertexPositions[he.vertex()];
    }
    /**
     * @brief Computes the volume variation vector for a given halfedge.
     * 
     * @param vpg The vertex position geometry containing face normals and areas.
     * @param he The halfedge for which to compute the volume variation vector.
     * @return Vector3 The volume variation vector.
     */
    
    inline Vector3 computeHalfedgeVolumeVariationVector(geometrycentral::surface::VertexPositionGeometry &vpg, GCHalfedge he)
    {
        std::size_t fID = he.face().getIndex();
        bool interiorHalfedge = he.isInterior();
        Vector3 volGrad{0, 0, 0};
        if (interiorHalfedge)
        {
            volGrad = vpg.faceNormals[fID] * vpg.faceAreas[fID] / 3;
        }
        return volGrad;
    }

    /**
     * @brief Computes the mean curvature vector for a given halfedge.
     * 
     * @param vpg The vertex position geometry containing face normals.
     * @param he The halfedge for which to compute the mean curvature vector.
     * @return Vector3 The mean curvature vector.
     */
    inline Vector3 computeHalfedgeMeanCurvatureVector(geometrycentral::surface::VertexPositionGeometry &vpg, GCHalfedge he)
    {
        std::size_t fID = he.face().getIndex();
        std::size_t fID_he_twin = he.twin().face().getIndex();
        bool interiorHalfedge = he.isInterior();
        bool interiorTwinHalfedge = he.twin().isInterior();
        Vector3 areaGrad{0, 0, 0};
        if (interiorHalfedge)
        {
            areaGrad += 0.25 * cross(vpg.faceNormals[fID], vecFromHalfedge(he.next(), vpg));
        }
        if (interiorTwinHalfedge)
        {
            areaGrad += 0.25 * cross(vpg.faceNormals[fID_he_twin], vecFromHalfedge(he.twin().next().next(), vpg));
        }
        return areaGrad / 2;
    }
    /**
     * @brief Computes the Gaussian curvature vector for a given halfedge.
     * 
     * @param vpg The vertex position geometry containing edge dihedral angles.
     * @param he The halfedge for which to compute the Gaussian curvature vector.
     * @return Vector3 The Gaussian curvature vector.
     */
    inline Vector3 computeHalfedgeGaussianCurvatureVector(geometrycentral::surface::VertexPositionGeometry &vpg, GCHalfedge he)
    {
        Vector3 gaussVec{0, 0, 0};
        if (!he.edge().isBoundary())
        {
            gaussVec = 0.5 * vpg.edgeDihedralAngles[he.edge()] * (-vecFromHalfedge(he, vpg)).unit();
        }
        return gaussVec;
    }
    /**
     * @brief Computes the gradient of the dihedral angle for a given halfedge and vertex.
     * 
     * @param geom The vertex position geometry containing edge lengths and face normals.
     * @param he The halfedge for which to compute the gradient.
     * @param v The vertex associated with the halfedge.
     * @return Vector3 The gradient of the dihedral angle.
     */
    inline Vector3 dihedralAngleGradient(geometrycentral::surface::VertexPositionGeometry &geom, GCHalfedge he, GCVertex v)
    {
        // std::cout <<"here arrived"<<std::endl;
        double l = geom.edgeLengths[he.edge()];
        if (he.edge().isBoundary())
        {
            return Vector3{0, 0, 0};
        }
        else if (he.vertex() == v)
        {
            return (geom.halfedgeCotanWeights[he.next().next()] *
                        geom.faceNormals[he.face()] +
                    geom.halfedgeCotanWeights[he.twin().next()] *
                        geom.faceNormals[he.twin().face()]) /
                   l;
        }
        else if (he.next().vertex() == v)
        {
            return (geom.halfedgeCotanWeights[he.twin().next().next()] *
                        geom.faceNormals[he.twin().face()] +
                    geom.halfedgeCotanWeights[he.next()] *
                        geom.faceNormals[he.face()]) /
                   l;
        }
        else if (he.next().next().vertex() == v)
        {
            return (-(geom.halfedgeCotanWeights[he.next().next()] +
                      geom.halfedgeCotanWeights[he.next()]) *
                    geom.faceNormals[he.face()]) /
                   l;
        }
        else
        {
            return Vector3{0, 0, 0};
        }
    }
    /**
     * @brief Computes the gradient of the area for a given vertex.
     * 
     * @param vpg The vertex position geometry containing face normals.
     * @param v The vertex for which to compute the area gradient.
     * @return Vector3 The area gradient at the vertex.
     */
    inline Vector3 computeGradAreaPointBased(geometrycentral::surface::VertexPositionGeometry &vpg, GCVertex &v)
    {
        Vector3 at_xm{0, 0, 0};

        for (GCHalfedge he : v.outgoingHalfedges())
        {
            assert(he.vertex() == v); // true
            GCFace f = he.face();
            // facenormals[index] may be faster, check the indexing speed next
            at_xm += 0.5 * cross(vpg.faceNormal(f),
                                 vecFromHalfedge(he.next(), vpg));
        }
        return at_xm; // SHOULD HAVE 1/3.0 but it does not work with FD computations
    }
    /**
     * @brief Computes the value of the Helfrich energy.
     * 
     * @return double The computed Helfrich energy value.
     */
    double HelfrichEnergyADE::Value()
    {
        geom->requireVertexDualAreas();
        geom->requireVertexMeanCurvatures();

        // parameters
        double kb = 1;
        //double alpha = 2 / geometrycentral::PI;
        double alpha = _alpha;
        double dela0 = targetValue;


        double m0 = 0;
        double m1 = 0;
        double m2 = 0;
        m0 = (geom->vertexDualAreas.raw().array()).sum();
        m1 = (geom->vertexMeanCurvatures.raw().array()).sum();
        m2 = (geom->vertexDualAreas.raw().array() * ((geom->vertexMeanCurvatures.raw().array() /
                                                      geom->vertexDualAreas.raw().array()) *
                                                     (geom->vertexMeanCurvatures.raw().array() /
                                                      geom->vertexDualAreas.raw().array())))
                 .sum();

        double radius = std::sqrt(m0/4/geometrycentral::PI);
        double dela0D = 4 * geometrycentral::PI * radius * dela0;
        double H0 = 0;         
        // std::cout << "  * Energy is calculated-based on new formulation" << std::endl;
        // std::cout << "  * Energy is" << m2*weight <<std::endl;
        double energy = 2 * kb * m2 + (2 * alpha * kb * geometrycentral::PI) * m1 * m1 / m0 -
                        4 * kb * H0 * m1 - 2 * geometrycentral::PI * alpha * kb * (dela0D)*m1 / m0 +
                        2 * kb * H0 * H0 * m0 + (0.5) * alpha * kb * geometrycentral::PI * (dela0D) * (dela0D) * (1 / m0);
        //double energy = 2 * kb * m2;
        return energy * weight;
    }
    /**
     * @brief Computes the differential of the Helfrich energy and stores it in the output matrix.
     * 
     * @param output The matrix to store the computed differential.
     */
    void HelfrichEnergyADE::Differential(Eigen::MatrixXd &output)
    {

        // std::cout << "  * Derivative" << std::endl;

        double kb = 1;
        //double alpha = 2 / geometrycentral::PI;
        double alpha = _alpha;
        double dela0 = targetValue;


        VertexIndices inds = mesh->getVertexIndices();

        geom->requireVertexDualAreas();
        geom->requireVertexMeanCurvatures();
        geom->requireEdgeLengths();
        geom->requireHalfedgeCotanWeights();
        geom->requireFaceNormals();
        double m0 = 0;
        double m1 = 0;
        double m2 = 0;
        m0 = (geom->vertexDualAreas.raw().array()).sum();
        // std::cout << m0<< std::endl;

        m1 = (geom->vertexMeanCurvatures.raw().array()).sum();
        m2 = (geom->vertexDualAreas.raw().array() * ((geom->vertexMeanCurvatures.raw().array() /
                                                      geom->vertexDualAreas.raw().array()) *
                                                     (geom->vertexMeanCurvatures.raw().array() /
                                                      geom->vertexDualAreas.raw().array())))
                 .sum();

        double radius = std::sqrt(m0/4/geometrycentral::PI);
        double dela0D = 4 * geometrycentral::PI * radius * dela0;
        double H0 = 0;

        for (std::size_t i = 0; i < mesh->nVertices(); ++i)
        {
            GCVertex v_i = mesh->vertex(i);
            Vector3 m0_xm = computeGradAreaPointBased(*geom, v_i); // Area Gradient
            Vector3 m1_xm{0, 0, 0};
            Vector3 m2_xm{0, 0, 0};
            Vector3 bendingForceVec{0, 0, 0};
            Vector3 bendingForceVec_areaGrad{0, 0, 0};
            Vector3 bendingForceVec_gaussVec{0, 0, 0};
            Vector3 bendingForceVec_schlafliVec{0, 0, 0};
            Vector3 vi_xm{0, 0, 0};
            Vector3 nu_xm{0, 0, 0};
            Vector3 dela_xm{0, 0, 0};

            // std::cout << "  mean curvature " <<geom->vertexMeanCurvatures.raw().array() <<std::endl;

            double Hi = geom->vertexMeanCurvatures[i] / geom->vertexDualAreas[i];
            double Ai = geom->vertexDualAreas[i];

            double H0i = 0.0;
            double Kbi = 1.0;
            double Kdi = 1.0;

            for (GCHalfedge he : v_i.outgoingHalfedges())
            {
                std::size_t fID = he.face().getIndex();
                // Initialize local variables for computation
                std::size_t i_vj = he.tipVertex().getIndex();
                double Hj = geom->vertexMeanCurvatures[i_vj] / geom->vertexDualAreas[i_vj];
                double Aj = geom->vertexDualAreas[i_vj];

                double H0j = 0.0;
                double Kbj = 1.0;
                double Kdj = 1.0;
                bool interiorHalfedge = he.isInterior();
                bool boundaryEdge = he.edge().isBoundary();
                bool boundaryNeighborVertex = he.next().vertex().isBoundary();

                Vector3 areaGrad = 2 * computeHalfedgeMeanCurvatureVector(*geom, he);

                Vector3 gaussVec = computeHalfedgeGaussianCurvatureVector(*geom, he);
                Vector3 schlafliVec1;
                Vector3 schlafliVec2;
                schlafliVec1 = geom->edgeLengths[he.edge()] * dihedralAngleGradient(*geom, he, he.vertex());
                schlafliVec2 = geom->edgeLengths[he.twin().edge()] *
                                   dihedralAngleGradient(*geom, he.twin(), he.vertex()) +
                               geom->edgeLengths[he.next().edge()] *
                                   dihedralAngleGradient(*geom, he.next(), he.vertex()) +
                               geom->edgeLengths[he.twin().next().next().edge()] *
                                   dihedralAngleGradient(*geom, he.twin().next().next(), he.vertex());

                bendingForceVec_schlafliVec -= (Kbi * (Hi - H0i) * schlafliVec1 + Kbj * (Hj - H0j) * schlafliVec2);
                bendingForceVec_areaGrad -= (Kbi * (H0i * H0i - Hi * Hi) / 3 +
                                             Kbj * (H0j * H0j - Hj * Hj) * 2 / 3) *
                                            areaGrad;
                bendingForceVec_gaussVec -=
                    (Kbi * (Hi - H0i) + Kbj * (Hj - H0j)) * gaussVec;

                m1_xm += gaussVec + 0.5 * (schlafliVec1 + schlafliVec2); // Calculated by me!

                m2_xm += (Hi * schlafliVec1 + Hj * schlafliVec2) +
                         ((-Hi * Hi) / 3 + (-Hj * Hj) * 2 / 3) * areaGrad +
                         (Hi + Hj) * gaussVec; // Check Mem3DG paper
            }

            auto Ex = 2 * kb * m2_xm + m1_xm * (alpha * kb * geometrycentral::PI * 4 * m1 / m0 - 4 * kb * H0 - 2 * geometrycentral::PI * alpha * kb * dela0D / m0) +
                      m0_xm * (2 * kb * H0 * H0 - 2 * alpha * kb * geometrycentral::PI * m1 * m1 / m0 / m0 + 2 * alpha * kb * geometrycentral::PI * dela0D * m1 / m0 / m0 - (0.5) * alpha * kb * geometrycentral::PI * dela0D * dela0D / m0 / m0);
            //auto Ex = 2 * kb * m2_xm;
            MatrixUtils::addToRow(output, inds[v_i], weight * Ex);
        }
    } // Differential


    /**
     * @brief Gets the exponents of the Helfrich energy.
     * 
     * @return Vector2 The exponents of the energy.
     */
    Vector2 HelfrichEnergyADE::GetExponents()
    {
        return Vector2{1, 0};
    }
    /**
     * @brief Adds the bi-Laplacian metric term to the provided list of metric terms.
     * 
     * @param terms The vector to which the metric term will be added.
     */
    void HelfrichEnergyADE::AddMetricTerm(std::vector<MetricTerm *> &terms)
    {
        std::cout << "  * Adding bi-Laplacian to metric for Helfrich energy" << std::endl;
        BiLaplacianMetricTerm *term = new BiLaplacianMetricTerm(mesh, geom);
        terms.push_back(term);
    }

    /**
     * @brief Resets the target values for the Helfrich energy.
     * 
     * @param newTargetValue The new target value for the reduced volume.
     * @param newWeight The new weight for the energy.
     * @param newalpha The new alpha parameter for the energy.
     */
    void HelfrichEnergyADE::ResetTargetsHelfrich(double newTargetValue, double newWeight, double newalpha)
    {
        // std::cout<<" Target value for the reduced volume is changed to: "<< newTargetValue<<std::endl;
        targetValue = newTargetValue;
        weight = newWeight;
        _alpha = newalpha;
    }

} // namespace biorsurfaces
