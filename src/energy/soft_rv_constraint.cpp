#include "energy/soft_rv_constraint.h"
#include "matrix_utils.h"
#include "surface_derivatives.h"

namespace biorsurfaces
{
    SoftReducedVolumeConstraint::SoftReducedVolumeConstraint(MeshPtr mesh_, GeomPtr geom_, double weight_, double targetValue_)
    {
        mesh = mesh_;
        geom = geom_;
        weight = weight_;
        targetValue = targetValue_;
        initialVolume = totalVolume(geom, mesh);
    }
inline Vector3 vecFromHalfedge(const  GCHalfedge  he, const geometrycentral::surface::VertexPositionGeometry &vpg) {
    return vpg.inputVertexPositions[he.next().vertex()] -
         vpg.inputVertexPositions[he.vertex()];
}

inline Vector3 computeHalfedgeVolumeVariationVector(geometrycentral::surface::VertexPositionGeometry &vpg, GCHalfedge he) {

  std::size_t fID = he.face().getIndex();
  bool interiorHalfedge = he.isInterior();
  Vector3 volGrad{0, 0, 0};
  if (interiorHalfedge) {
    volGrad = vpg.faceNormals[fID] * vpg.faceAreas[fID] / 3;
  }
  return volGrad;
}

inline Vector3  computeHalfedgeMeanCurvatureVector(geometrycentral::surface::VertexPositionGeometry &vpg, GCHalfedge he) {
        std::size_t fID = he.face().getIndex();
        std::size_t fID_he_twin = he.twin().face().getIndex();
        bool interiorHalfedge = he.isInterior();
        bool interiorTwinHalfedge = he.twin().isInterior();
        Vector3 areaGrad{0, 0, 0};
        if (interiorHalfedge) {
            areaGrad += 0.25 * cross(vpg.faceNormals[fID], vecFromHalfedge(he.next(),vpg));
        }
        if (interiorTwinHalfedge) {
            areaGrad += 0.25 * cross(vpg.faceNormals[fID_he_twin], vecFromHalfedge(he.twin().next().next(),vpg));
        }
        return areaGrad / 2;
    }

    inline Vector3 computeHalfedgeGaussianCurvatureVector(geometrycentral::surface::VertexPositionGeometry &vpg, GCHalfedge he) {
        Vector3 gaussVec{0, 0, 0};
        if (!he.edge().isBoundary()) {
            gaussVec = 0.5 * vpg.edgeDihedralAngles[he.edge()] *(-vecFromHalfedge(he,vpg)).unit();
        }
        return gaussVec;
}

    inline Vector3 dihedralAngleGradient(geometrycentral::surface::VertexPositionGeometry &geom,GCHalfedge he, GCVertex v) {
          //std::cout <<"here arrived"<<std::endl;
           double l = geom.edgeLengths[he.edge()];
           if (he.edge().isBoundary()) {
              return Vector3{0, 0, 0};
           } else if (he.vertex() == v) {
              return (   geom.halfedgeCotanWeights[he.next().next()] *
                         geom.faceNormals[he.face()] +
                         geom.halfedgeCotanWeights[he.twin().next()] *
                         geom.faceNormals[he.twin().face()]) /
                         l;
           } else if (he.next().vertex() == v) {
                 return (geom.halfedgeCotanWeights[he.twin().next().next()] *
                         geom.faceNormals[he.twin().face()] +
                         geom.halfedgeCotanWeights[he.next()] *
                         geom.faceNormals[he.face()]) /
                         l;
  } else if (he.next().next().vertex() == v) {
            return (-(geom.halfedgeCotanWeights[he.next().next()] +
                       geom.halfedgeCotanWeights[he.next()]) *
                       geom.faceNormals[he.face()]) /
                       l;
  } else {
            return Vector3{0, 0, 0};
  }
}

inline Vector3 computeGradAreaPointBased(geometrycentral::surface::VertexPositionGeometry &vpg, GCVertex &v){
     Vector3 at_xm{0,0,0};                                

for(GCHalfedge he : v.outgoingHalfedges()) {
  assert(he.vertex() == v); // true
  GCFace f = he.face();
  //facenormals[index] may be faster, check the indexing speed next
  at_xm += 0.5*cross(vpg.faceNormal(f),
                     vecFromHalfedge(he.next(),vpg));
  }
  return at_xm;//SHOULD HAVE 1/3.0 but it does not work with FD computations
}

inline double volumeDeviation(MeshPtr mesh, GeomPtr geom, double initialValue)
    {
        return (totalVolume(geom, mesh) - initialValue) / initialValue;
    }

    // Returns the current value of the energy.
    double SoftReducedVolumeConstraint::Value()
    {
        geom->requireVertexDualAreas();
        geom->requireVertexMeanCurvatures();
        double volume = totalVolume(geom, mesh);
        double m0 = 0;
        double m1 = 0;
        double m2 = 0; 
        m0 = (geom->vertexDualAreas.raw().array()).sum();
        m1 = (geom->vertexMeanCurvatures.raw().array()).sum();
        m2 = (geom->vertexDualAreas.raw().array()*((geom->vertexMeanCurvatures.raw().array() /
              geom->vertexDualAreas.raw().array())*(geom->vertexMeanCurvatures.raw().array() /
              geom->vertexDualAreas.raw().array()))).sum();
        return weight*(6*sqrt(geometrycentral::PI/m0)*(volume/m0)-targetValue)*(6*sqrt(geometrycentral::PI/m0)*(volume/m0)-targetValue);
    }

    // Returns the current differential of the energy, stored in the given
    // V x 3 matrix, where each row holds the differential (a 3-vector) with
    // respect to the corresponding vertex.
    void SoftReducedVolumeConstraint::Differential(Eigen::MatrixXd &output)
    {      
     // std::cout << "  * Derivative" << std::endl;
      VertexIndices inds = mesh->getVertexIndices();
      double volume = totalVolume(geom, mesh);

      geom->requireVertexDualAreas();
      geom->requireVertexMeanCurvatures();
      geom->requireEdgeLengths();
      geom->requireHalfedgeCotanWeights();
      geom->requireFaceNormals();

      for (std::size_t i = 0; i < mesh->nVertices(); ++i) {          
         GCVertex v_i = mesh->vertex(i);
         double m0 = 0;
         double m1 = 0;
         double m2 = 0; 

         m0 = (geom->vertexDualAreas.raw().array()).sum();
        //std::cout << m0<< std::endl;

         m1 = (geom->vertexMeanCurvatures.raw().array()).sum();
         m2 =(geom->vertexDualAreas.raw().array()*((geom->vertexMeanCurvatures.raw().array() /
              geom->vertexDualAreas.raw().array())*(geom->vertexMeanCurvatures.raw().array() /
              geom->vertexDualAreas.raw().array()))).sum();

         Vector3 m0_xm = computeGradAreaPointBased(*geom, v_i);//Area Gradient

         Vector3 m1_xm{0, 0, 0};
         Vector3 m2_xm{0, 0, 0};
         Vector3 bendingForceVec{0, 0, 0};
         Vector3 bendingForceVec_areaGrad{0, 0, 0};
         Vector3 bendingForceVec_gaussVec{0, 0, 0};
         Vector3 bendingForceVec_schlafliVec{0, 0, 0};
         Vector3 vi_xm{0, 0, 0};
         Vector3 nu_xm{0, 0, 0};
         Vector3 dela_xm{0, 0, 0};

        //std::cout << "  mean curvature " <<geom->vertexMeanCurvatures.raw().array() <<std::endl;

        double Hi = geom->vertexMeanCurvatures[i] / geom->vertexDualAreas[i];
        double Ai = geom->vertexDualAreas[i];

        double H0i = 0.0;
        double Kbi = 1.0;
        double Kdi = 1.0;

        for (GCHalfedge he : v_i.outgoingHalfedges()) {
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

            m1_xm +=   gaussVec+ 0.5*(schlafliVec1+schlafliVec2); //Calculated by me!
            m2_xm +=   (Hi * schlafliVec1 + Hj* schlafliVec2)+
                       ((-Hi * Hi) / 3 + (- Hj * Hj) * 2 / 3) * areaGrad+
                       (Hi + Hj) * gaussVec; // Check Mem3DG paper

            vi_xm += computeHalfedgeVolumeVariationVector(*geom, he);   
        }


        nu_xm = 2*6*sqrt(geometrycentral::PI/m0/m0/m0)*(vi_xm-(3/2.0)*m0_xm*(volume/m0))*(6*sqrt(geometrycentral::PI/m0)*(volume/m0)-targetValue);
        MatrixUtils::addToRow(output, inds[v_i], weight * nu_xm);
        }
    }

    // Get the exponents of this energy; only applies to tangent-point energies.
    Vector2 SoftReducedVolumeConstraint::GetExponents()
    {
        return Vector2{1, 0};
    }

    // Get a pointer to the current BVH for this energy.
    // Return 0 if the energy doesn't use a BVH.
    OptimizedClusterTree *SoftReducedVolumeConstraint::GetBVH()
    {
        return 0;
    }

    // Return the separation parameter for this energy.
    // Return 0 if this energy doesn't do hierarchical approximation.
    double SoftReducedVolumeConstraint::GetTheta()
    {
        return 0;
    }

    void SoftReducedVolumeConstraint::ResetTargetsRV(double newTargetValue,double newWeight) {
        //std::cout<<" Target value for the reduced volume is changed to: "<< newTargetValue<<std::endl;               
            targetValue = newTargetValue;
            weight = newWeight;
    }


} // namespace biorsurfaces