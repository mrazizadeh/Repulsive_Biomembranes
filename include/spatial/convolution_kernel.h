#pragma once

#include "rsurface_types.h"

namespace biorsurfaces
{
    class RieszKernel
    {
    public:
        RieszKernel(double s_);
        double Coefficient(Vector3 x, Vector3 y);

    private:
        double s;
    };
} // namespace biorsurfaces