#pragma once

#include "rsurface_types.h"

namespace biorsurfaces
{
    void writeMeshToOBJ(MeshPtr mesh, GeomPtr geom, GeomPtr geomOrig, bool writeAreaRatios, std::string output);
}
