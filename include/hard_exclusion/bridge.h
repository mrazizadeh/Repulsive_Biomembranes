#pragma once
#include "hard_exclusion/guard.h"
#include "rsurface_types.h"
#include "scene_file.h"

namespace biorsurfaces { namespace HardExclusion {
void ValidateScene(const scene::SceneData& data);
void Initialize(const scene::SceneData& data, const MeshPtr& mesh, const GeomPtr& geom);
bool Enabled();
double ConstraintRelativeError(const MeshPtr& mesh, const GeomPtr& geom);
double ConstraintRelativeTolerance();
hard_exclusion::Positions Capture(const MeshPtr& mesh, const GeomPtr& geom);
void Restore(const MeshPtr& mesh, const GeomPtr& geom, const hard_exclusion::Positions& x);
double Limit(const hard_exclusion::Positions& before, const hard_exclusion::Positions& proposed, const char* phase);
void RecordLinearSolve(int info, long iterations, double residual, double tolerance, bool finiteSolution, const char* stage);
void BeginStep(int step, const MeshPtr& mesh, const GeomPtr& geom, double objective);
void AcceptStep(const MeshPtr& mesh, const GeomPtr& geom, double objective);
void FailStep(const MeshPtr& mesh, const GeomPtr& geom, const std::string& reason);
void Finish();
}}
