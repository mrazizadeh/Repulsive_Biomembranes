#pragma once

#include "rsurface_types.h"
#include "matrix_utils.h"

namespace biorsurfaces
{
    namespace H2
    {
        void getTriplets(std::vector<Triplet> &triplets, const MeshPtr &mesh, const GeomPtr &geom, double epsilon);
    } // namespace H2
} // namespace biorsurfaces