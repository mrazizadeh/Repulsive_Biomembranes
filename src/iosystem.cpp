#include <cassert>
#include <iostream>
#include <netcdf>
#include <vector>
#include <Eigen/Core>

#include "iosystem.h"

namespace biorsurfaces{
    namespace nc = ::netCDF;
    using NcException = nc::exceptions::NcException;
    using NcFile = nc::NcFile;


    void TrajFile::writeCoords(
        const std::size_t idx,
        const Eigen::Matrix<double, Eigen::Dynamic, SPATIAL_DIMS, Eigen::RowMajor>
            &data) {
    if (!writeable)
        throw std::runtime_error("Cannot open an already opened ...");

    assert(data.rows() == nvertices_dim.getSize());

    coord_var.putVar({idx, 0, 0}, {1, nvertices_dim.getSize(), SPATIAL_DIMS},
                    data.data());
    }

    void TrajFile::writeTopoFrame(const std::size_t idx,
                                const Eigen::Matrix<std::uint32_t, Eigen::Dynamic,
                                                    3, Eigen::RowMajor> &data) {
    if (!writeable)
        throw std::runtime_error("Cannot open an already opened ...");

    assert(data.rows() == npolygons_dim.getSize());

    topo_frame_var.putVar(
        {idx, 0, 0}, {1, npolygons_dim.getSize(), POLYGON_ORDER}, data.data());
    }

    EigenVectorX3dr TrajFile::getCoords(const std::size_t idx) const {
    assert(idx < nFrames());

    EigenVectorX3dr vec(nvertices_dim.getSize(), SPATIAL_DIMS);
    coord_var.getVar({idx, 0, 0}, {1, nvertices_dim.getSize(), SPATIAL_DIMS},
                    vec.data());
    return vec;
    }

// topology frame
Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>
TrajFile::getTopoFrame(const std::size_t idx) const {
  Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor> vec(
      npolygons_dim.getSize(), POLYGON_ORDER);
  topo_frame_var.getVar(
      {idx, 0, 0}, {1, npolygons_dim.getSize(), POLYGON_ORDER}, vec.data());
  return vec;
}

// topology
Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor>
TrajFile::getTopology() const {
  Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 3, Eigen::RowMajor> vec(
      npolygons_dim.getSize(), POLYGON_ORDER);
  topology.getVar({0, 0}, {npolygons_dim.getSize(), POLYGON_ORDER}, vec.data());
  return vec;
}



}