#pragma once
#include <Eigen/Core>
#include <cmath>

namespace biorsurfaces { namespace HardExclusion {
inline bool LinearSolveAcceptable(int info,double residual,double tolerance,bool finiteSolution) {
    return info==int(Eigen::Success) && finiteSolution && std::isfinite(residual) && residual>=0 &&
           std::isfinite(tolerance) && tolerance>0 && residual<=tolerance;
}
}}
