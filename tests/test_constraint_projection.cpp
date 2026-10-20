#include "hard_exclusion/constraint_projection.h"
#include <iostream>
#include <string>

using namespace biorsurfaces::HardExclusion;
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
template <typename F>
void expectFailure(F f, const char* reason) {
    try { f(); } catch (const std::runtime_error& e) {
        require(std::string(e.what()) == reason, "unexpected failure reason");
        return;
    }
    throw std::runtime_error("failure was not propagated");
}
int main() {
    // A nonlinear volume constraint V=x^3. One Newton correction from x=1.1
    // leaves ~2.7% volume error, so a one-correction implementation must fail
    // this regression. Use a frozen Jacobian like the production Schur loop.
    double x = 1.1;
    const double jacobian = 3*x*x;
    int calls = 0;
    const int count = CorrectConstraintsToTolerance(
        [&]() { return RelativeAreaVolumeError(1, x*x*x, 1, 1); },
        [&](int) { x -= (x*x*x-1)/jacobian; ++calls; }, 0.001);
    require(count == calls && count > 1 && count <= MaxConstraintCorrections,
            "nonlinear correction did not iterate");
    require(std::abs(x*x*x-1) <= 0.001, "volume still violates tolerance");
    std::cout << "PASS nonlinear volume requires " << count
              << " corrections; relative residual=" << std::abs(x*x*x-1) << '\n';

    // Relative A/V checks must be scale-independent and enforce both targets.
    for (double scale : {0.001, 1., 1000.}) {
        const double a=10.026827372651271*scale*scale;
        const double v=1.0449299195772292*scale*scale*scale;
        require(RelativeAreaVolumeError(a, v*1.0030715, a, v)>0.001,
                "reported smoke-failure magnitude accepted");
        require(RelativeAreaVolumeError(a*1.0011, v, a, v)>0.001,
                "area-only failure accepted");
        require(RelativeAreaVolumeError(a*1.0005, v*0.9995, a, v)<=0.001,
                "admissible constraints rejected");
    }
    std::cout << "PASS separate area/volume checks and scale independence\n";

    int attempts=0;
    expectFailure([&]() {
        CorrectConstraintsToTolerance([]() { return 0.003; },
                                      [&](int) { ++attempts; },0.001);
    },"constraint_correction_budget_exhausted");
    require(attempts==MaxConstraintCorrections,"correction budget not enforced");
    std::cout << "PASS stalled correction stops at budget without accepting\n";

    attempts=0;
    expectFailure([&]() {
        CorrectConstraintsToTolerance(
            []() { return RelativeAreaVolumeError(1,-1,1,1); },
            [&](int) { ++attempts; },0.001);
    },"invalid_area_or_volume_during_correction");
    require(attempts==0,"invalid geometry sent to correction");
    expectFailure([&]() {
        CorrectConstraintsToTolerance([]() { return 0.01; },
            [](int) { throw std::runtime_error("collision_rejection"); },0.001);
    },"collision_rejection");
    std::cout << "PASS invalid geometry and collision rejection propagate\n";

    attempts=0;
    const int zero=CorrectConstraintsToTolerance([]() { return 0.; },
                                                [&](int) { ++attempts; },0.001);
    require(zero==1 && attempts==1,"missing checked zero-correction phase");
    std::cout << "PASS initially satisfied constraints still check correction phase\n";
}
