#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace biorsurfaces { namespace HardExclusion {
constexpr int MaxConstraintCorrections = 20;

inline double RelativeAreaVolumeError(double area, double volume,
                                      double targetArea, double targetVolume) {
    if (!std::isfinite(area) || !std::isfinite(volume) ||
        !std::isfinite(targetArea) || !std::isfinite(targetVolume) ||
        area <= 0 || volume <= 0 || targetArea <= 0 || targetVolume <= 0)
        return std::numeric_limits<double>::infinity();
    return std::max(std::abs(area / targetArea - 1),
                    std::abs(volume / targetVolume - 1));
}

// The correction callback must guard each motion and refresh geometry before
// error() is called again. Keep at least one checked correction phase, including
// a zero correction, for the complete-step phase audit. The error is measured
// from geometry, not the mixed-unit Schur right-hand side.
template <typename Error, typename Correction>
int CorrectConstraintsToTolerance(Error error, Correction correct,
                                 double tolerance,
                                 int limit = MaxConstraintCorrections) {
    if (!std::isfinite(tolerance) || tolerance <= 0 || limit < 1)
        throw std::runtime_error("invalid_constraint_correction_policy");
    for (int n = 0; ; ++n) {
        const double residual = error();
        if (!std::isfinite(residual) || residual < 0)
            throw std::runtime_error("invalid_area_or_volume_during_correction");
        if (n > 0 && residual <= tolerance) return n;
        if (n == limit)
            throw std::runtime_error("constraint_correction_budget_exhausted");
        correct(n);
    }
}
}}
