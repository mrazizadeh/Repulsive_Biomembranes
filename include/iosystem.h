#pragma once

#include <exception>
#include <iostream>
#include <netcdf>
#include <vector>
#include "trajfile_constants.h"

#include <Eigen/Core>

#include <geometrycentral/surface/halfedge_mesh.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>

namespace biorsurfaces {

namespace gc = ::geometrycentral;
namespace gcs = ::geometrycentral::surface;

namespace nc = ::netCDF;
using EigenVectorX3dr =   Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;

class TrajFile {

private:

  TrajFile *operator=(const TrajFile &rhs) = delete;
  TrajFile(const TrajFile &rhs) = delete;

public:
  using NcFile = nc::NcFile;
  using NcException = nc::exceptions::NcException;

  TrajFile() : fd(nullptr), writeable(false){};

  void open(const std::string &filename, const NcFile::FileMode fMode) {
    if ((fd != nullptr) && (fMode != NcFile::read)) {
        throw std::runtime_error("Cannot open an already opened ...");
    }
    fd = new NcFile(filename, fMode);
    writeable = fMode != NcFile::read;

    frame_dim = fd->getDim(FRAME_NAME);
    npolygons_dim = fd->getDim(NPOLYGONS_NAME);
    nvertices_dim = fd->getDim(NVERTICES_NAME);
    ncorners_dim = fd->getDim(NCORNERS_NAME);
    spatial_dim = fd->getDim(SPATIAL_DIMS_NAME);
    polygon_order_dim = fd->getDim(POLYGON_ORDER_NAME);
    topology = fd->getVar(TOPO_VAR);
  }

  void createNewFile(const std::string &filename, gcs::SurfaceMesh &mesh,
                     gcs::VertexPositionGeometry &refVpg,
                     const NcFile::FileMode fMode) {
    if (fd != nullptr) {
        throw std::runtime_error("Cannot open an already opened ...");
    }

    writeable = true;

    fd = new NcFile(filename, fMode);
    // initialize data
    fd->putAtt(CONVENTIONS_NAME, CONVENTIONS_VALUE);
    fd->putAtt(CONVENTIONS_VERSION_NAME, CONVENTIONS_VERSION_VALUE);

    frame_dim = fd->addDim(FRAME_NAME);
    npolygons_dim = fd->addDim(NPOLYGONS_NAME, mesh.nFaces());
    nvertices_dim = fd->addDim(NVERTICES_NAME, mesh.nVertices());
    ncorners_dim = fd->addDim(NCORNERS_NAME, mesh.nCorners());
    spatial_dim = fd->addDim(SPATIAL_DIMS_NAME, SPATIAL_DIMS);
    polygon_order_dim = fd->addDim(POLYGON_ORDER_NAME, POLYGON_ORDER);

    // Initialize topology data block
    topology = fd->addVar(TOPO_VAR, nc::ncUint, {npolygons_dim, polygon_order_dim});

    // Populate topology data
    Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor> faceMatrix{
        mesh.getFaceVertexMatrix<std::uint32_t>()};

    std::uint32_t *topodata = faceMatrix.data();
    topology.putVar(topodata);

    auto mapped_position_matrix = gc::EigenMap<double, 3, Eigen::RowMajor>(refVpg.inputVertexPositions);

    // Ensure that mapped_position_matrix contains the correct data
    std::cout << "Mapped Position Matrix: " << std::endl << mapped_position_matrix.data() << std::endl;
    // Get a pointer to the data
    double *coorddata = mapped_position_matrix.data();
    // Add the variable to the NetCDF file
    coord_var = fd->addVar(COORD_VAR, netCDF::ncDouble, {frame_dim, nvertices_dim, spatial_dim});

    coord_var.putVar(coorddata);

    topo_frame_var = fd->addVar(TOPO_FRAME_VAR, netCDF::ncUint,
                                {frame_dim, npolygons_dim, polygon_order_dim});
  }

  TrajFile(TrajFile &&rhs) = default;

  static TrajFile newFile(const std::string &filename, gcs::SurfaceMesh &mesh,
                          gcs::VertexPositionGeometry &refVpg,
                          bool replace = false) {
    if (replace)
      return TrajFile(filename, mesh, refVpg, NcFile::replace);
    else
      return TrajFile(filename, mesh, refVpg, NcFile::newFile);
  };

  static TrajFile openRW(const std::string &filename) {
    return TrajFile(filename, NcFile::write);
  }

  static TrajFile openReadOnly(const std::string &filename) {
    return TrajFile(filename, NcFile::read);
  };


  void getNcFrame(int &frame) const {
    int maxFrame = nFrames() - 1;
    if (frame > maxFrame || frame < -(maxFrame + 1)) {
        throw std::runtime_error("Snapshot frame exceed limiting frame index!");
    } else if (frame < 0) {
      frame = frame + maxFrame + 1;
    }
  }

  /**
   * @brief Destructor frees the bound NcFile
   */
  ~TrajFile() { delete fd; };

  /**
   * @brief Check if the bound NcFile was opened in write mode
   *
   * @return True if writable
   */
  inline bool isWriteable() { return writeable; };

  inline std::size_t nFrames() const { return frame_dim.getSize(); };

  /**
   * @brief Validate whether or not the metadata follows convention
   *
   * @return True if okay. False otherwise.
   */



  void writeCoords(const std::size_t idx, const EigenVectorX3dr &data);

  void writeTopoFrame(const std::size_t idx,
                      const Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3,
                                          Eigen::RowMajor> &data);



  EigenVectorX3dr getCoords(const std::size_t idx) const;

  Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>
  getTopoFrame(const std::size_t idx) const;

  Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>
  getTopology() const;

private:
  TrajFile(const std::string &filename, const NcFile::FileMode fMode)
      : filename(filename), // fd(new NcFile(filename, fMode)),
        writeable(fMode != NcFile::read) {
    open(filename, fMode);
  }

  TrajFile(const std::string &filename, gcs::SurfaceMesh &mesh,
           gcs::VertexPositionGeometry &refVpg, const NcFile::FileMode fMode)
      : filename(filename), // fd(new NcFile(filename, fMode)),
        writeable(true) {
    createNewFile(filename, mesh, refVpg, fMode);
  }

  /// Bound NcFile
  NcFile *fd;

  // Save dimensions
  nc::NcDim frame_dim;
  nc::NcDim npolygons_dim;
  nc::NcDim nvertices_dim;
  nc::NcDim ncorners_dim;
  nc::NcDim spatial_dim;
  nc::NcDim polygon_order_dim;
  nc::NcVar topo_frame_var;
  // Save variables
  nc::NcVar topology;
  nc::NcVar coord_var;
  /// Filepath to file
  std::string filename;
  /// Writeable status
  bool writeable = true;
};
} // namespace solver
