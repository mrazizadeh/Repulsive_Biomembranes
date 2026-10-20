#include "hard_exclusion/bridge.h"
#include "hard_exclusion/linear_solve.h"
#include "hard_exclusion/constraint_projection.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace biorsurfaces { namespace HardExclusion {
namespace {
std::unique_ptr<hard_exclusion::Guard> guard;
hard_exclusion::Positions stepStart;
std::ofstream motions, steps, linearSolves;
double targetArea=0, targetVolume=0, relativeTolerance=0, initialEdge=0;
double startEnergy=0, smallestFraction=1;
int currentStep=0, acceptedSteps=0, motionCalls=0, limitedCalls=0, phaseMask=0;
std::string jsonString(const std::string& s) {
    std::string out="\"";
    for (char c:s) { if(c=='\\' || c=='\"') out+='\\'; if(c=='\n' || c=='\r') out+=' '; else out+=c; }
    return out+'\"';
}
void status(const char* state, const std::string& reason) {
    std::ofstream out("hard_exclusion_status.json");
    out << "{\"status\":" << jsonString(state) << ",\"reason\":" << jsonString(reason)
        << ",\"accepted_steps\":" << acceptedSteps << ",\"attempted_step\":" << currentStep << "}\n";
}
void require(bool condition,const char* message) { if(!condition) throw std::runtime_error(message); }
void writeDiagnosticMesh(const char* name,const hard_exclusion::Positions& x,const hard_exclusion::Faces& f) {
    std::ofstream out(name);
    if(!out) { std::cerr << "Cannot write diagnostic mesh: " << name << std::endl; return; }
    out << std::setprecision(17);
    for(int i=0;i<x.rows();++i) out << "v " << x(i,0) << ' ' << x(i,1) << ' ' << x(i,2) << '\n';
    for(int i=0;i<f.rows();++i) out << "f " << f(i,0)+1 << ' ' << f(i,1)+1 << ' ' << f(i,2)+1 << '\n';
}
void checkConstraints(const hard_exclusion::Positions& x) {
    const double a=hard_exclusion::area(x,guard->faces()), v=hard_exclusion::volume(x,guard->faces());
    require(std::isfinite(a) && std::isfinite(v) && a>0 && v>0,"invalid_area_or_volume");
    require(RelativeAreaVolumeError(a,v,targetArea,targetVolume)<=relativeTolerance,
            "area_or_volume_outside_declared_relative_tolerance");
}
}
bool Enabled() { return bool(guard); }
double ConstraintRelativeError(const MeshPtr& mesh,const GeomPtr& geom) {
    require(Enabled(),"constraint_error_requested_without_hard_exclusion");
    const auto x=Capture(mesh,geom);
    return RelativeAreaVolumeError(hard_exclusion::area(x,guard->faces()),
                                   hard_exclusion::volume(x,guard->faces()),targetArea,targetVolume);
}
double ConstraintRelativeTolerance() { return relativeTolerance; }
hard_exclusion::Positions Capture(const MeshPtr& mesh,const GeomPtr& geom) {
    hard_exclusion::Positions x(mesh->nVertices(),3);
    auto ids=mesh->getVertexIndices();
    for(GCVertex v:mesh->vertices()) { Vector3 p=geom->inputVertexPositions[v]; x.row(ids[v]) << p.x,p.y,p.z; }
    return x;
}
void Restore(const MeshPtr& mesh,const GeomPtr& geom,const hard_exclusion::Positions& x) {
    require(x.rows()==static_cast<int>(mesh->nVertices()),"hard_exclusion_connectivity_changed");
    auto ids=mesh->getVertexIndices();
    for(GCVertex v:mesh->vertices()) { int i=ids[v]; geom->inputVertexPositions[v]=Vector3{x(i,0),x(i,1),x(i,2)}; }
    geom->refreshQuantities();
}
void ValidateScene(const scene::SceneData& d) {
    if(!d.hardExclusion) return;
    require(!d.remeshEnabled,"Hard exclusion requires remesh off in BOTH matched cases");
    require(d.defaultMethod==GradientMethod::HsProjectedIterative,"Hard exclusion currently supports method hs only");
    require(!d.gradual_rv && !d.autoComputeVolumeTarget && !d.allowBarycenterShift && !d.disableNearField,
            "Hard exclusion does not support ramps or altered metric/constraint settings");
    require(d.obstacles.empty() && d.implicitBarriers.empty() && d.vertexPins.empty() && d.vertexNormals.empty(),
            "Hard exclusion pilot supports only self exclusion of one closed membrane");
    require(d.potentials.size()==1 && d.potentials[0].type==scene::PotentialType::HelfrichADE && d.potentials[0].weight==1,
            "Hard exclusion pilot needs exactly minimize helfrichADE 1");
    require(d.constraints.size()==3 && d.constraints[0].type==scene::ConstraintType::Barycenter &&
            d.constraints[1].type==scene::ConstraintType::TotalVolume && d.constraints[2].type==scene::ConstraintType::ReducedVolume,
            "Use constraints in order: barycenter then volume 1 then reduced_volume 1");
    for(const auto& c:d.constraints) require(c.numIterations==0 && c.targetMultiplier==1 && c.targetAddition==0,
            "Hard exclusion needs constant unscaled constraint targets");
    require(std::isfinite(d.reducedVolume) && d.reducedVolume>0 && d.reducedVolume<=1 &&
            std::isfinite(d.weightHelfrich) && d.weightHelfrich>0 && std::isfinite(d.alphaADE) && d.alphaADE>=0 &&
            std::isfinite(d.reducedArea) && d.netCDFWriteIndex==1 && d.iterationLimit>0,
            "Invalid model parameters or write_index not equal to 1");
    require(std::isfinite(d.hardMinRelative) && d.hardMinRelative>0 && std::isfinite(d.hardToleranceRelative) &&
            d.hardToleranceRelative>0 && d.hardToleranceRelative<d.hardMinRelative &&
            std::isfinite(d.hardConstraintTolerance) && d.hardConstraintTolerance>0 && d.hardConstraintTolerance<=0.001,
            "Invalid hard-exclusion separation/tolerance or A/V tolerance above 0.001");
}
void Initialize(const scene::SceneData& d,const MeshPtr& mesh,const GeomPtr& geom) {
    if(!d.hardExclusion) return;
    const char* checked=std::getenv("HARD_EXCLUSION_INPUT_AUDITED");
    require(checked && std::string(checked)=="1",
            "Use audited runner: static VF/EE checks alone cannot validate an arbitrary initially intersecting mesh");
    ValidateScene(d);
    auto x=Capture(mesh,geom);
    hard_exclusion::Faces f=mesh->getFaceVertexMatrix<int>();
    hard_exclusion::validate_closed_manifold(x,f);
    initialEdge=hard_exclusion::mean_edge_length(x,f);
    hard_exclusion::Options o;
    o.minimum_separation=d.hardMinRelative*initialEdge; o.tolerance=d.hardToleranceRelative*initialEdge;
    guard.reset(new hard_exclusion::Guard(f,o));
    auto r=guard->check(x,x);
    if(!r.safe) { status("invalid_input",r.reason); throw std::runtime_error("Hard exclusion: initially invalid geometry: "+r.reason); }
    targetVolume=hard_exclusion::volume(x,f);
    targetArea=std::pow(6*std::sqrt(geometrycentral::PI)*targetVolume/d.reducedVolume,2./3.);
    relativeTolerance=d.hardConstraintTolerance;
    try { checkConstraints(x); } catch(const std::exception& e) { status("invalid_input",e.what()); throw; }
    std::ofstream config("hard_exclusion_config.json"); config << std::setprecision(17)
        << "{\"method\":\"Hs_fixed_connectivity_TightInclusion_BFS\",\"mean_initial_edge\":" << initialEdge
        << ",\"minimum_separation\":" << o.minimum_separation << ",\"ccd_tolerance\":" << o.tolerance
        << ",\"target_area\":" << targetArea << ",\"target_volume\":" << targetVolume
        << ",\"constraint_relative_tolerance\":" << relativeTolerance
        << ",\"constraint_correction_max_iterations\":" << MaxConstraintCorrections
        << ",\"constraint_correction_method\":\"guarded_frozen_Schur_until_relative_AV_tolerance\""
        << ",\"objective_increase_relative_tolerance\":1e-8,\"repulsion\":\"" << (d.repulsionOff?"off":"on")
        << "\",\"remesh\":false,\"ccd_commit\":\"6f84001a790a9d3362b7c13d410037140d83fa2c\"}\n";
    motions.open("hard_exclusion_motions.csv"); steps.open("hard_exclusion_steps.csv");
    linearSolves.open("hard_exclusion_linear_solves.csv");
    require(bool(motions) && bool(steps) && bool(linearSolves) && bool(config),"cannot_open_hard_exclusion_audit_files");
    linearSolves << "step,stage,eigen_info,iterations,reported_residual,tolerance,finite_solution,accepted\n" << std::setprecision(17);
    motions << "step,phase,safe,fraction,backtracks,face_pairs,vf_queries,ee_queries,reason\n" << std::setprecision(17);
    steps << "step,area,volume,area_relative_error,volume_relative_error,objective,objective_change,rms_displacement,motion_calls,limited_calls,minimum_fraction\n" << std::setprecision(17);
    status("initialized","input_validated");
    std::cout << "*** HARD EXCLUSION ACTIVE: continuous descent and correction checks; fixed connectivity ***" << std::endl;
}
double Limit(const hard_exclusion::Positions& before,const hard_exclusion::Positions& proposed,const char* phase) {
    if(!Enabled()) return 1;
    auto r=guard->filter(before,proposed); ++motionCalls;
    if(std::string(phase)=="descent") phaseMask|=1;
    if(std::string(phase)=="constraint_correction") phaseMask|=2;
    if(std::string(phase)=="barycenter_projection") phaseMask|=4;
    if(r.fraction<1) ++limitedCalls;
    smallestFraction=std::min(smallestFraction,r.fraction);
    motions << currentStep << ',' << phase << ',' << r.safe << ',' << r.fraction << ',' << r.backtracks << ','
            << r.face_pairs << ',' << r.vertex_face_queries << ',' << r.edge_edge_queries << ',' << r.reason << std::endl;
    require(r.safe && r.fraction>0,"hard_exclusion_no_admissible_motion");
    return r.fraction;
}
void RecordLinearSolve(int info,long iterations,double residual,double tolerance,bool finiteSolution,const char* stage) {
    if(!Enabled()) return;
    const bool accepted=LinearSolveAcceptable(info,residual,tolerance,finiteSolution);
    linearSolves << currentStep << ',' << stage << ',' << info << ',' << iterations << ',' << residual << ','
                 << tolerance << ',' << finiteSolution << ',' << accepted << std::endl;
    if(!accepted) {
        std::ostringstream message;
        message << std::setprecision(17) << "gmres_failed:stage=" << stage << ":info=" << info << ":iterations=" << iterations
                << ":reported_residual=" << residual << ":tolerance=" << tolerance;
        throw std::runtime_error(message.str());
    }
}
void BeginStep(int step,const MeshPtr& mesh,const GeomPtr& geom,double objective) {
    if(!Enabled()) return;
    currentStep=step; stepStart=Capture(mesh,geom); startEnergy=objective;
    motionCalls=limitedCalls=phaseMask=0; smallestFraction=1;
    require(std::isfinite(objective),"nonfinite_start_objective");
    status("running","step_in_progress");
}
void AcceptStep(const MeshPtr& mesh,const GeomPtr& geom,double objective) {
    if(!Enabled()) return;
    require(phaseMask==7,"missing_guarded_phase_in_complete_step");
    auto f=mesh->getFaceVertexMatrix<int>();
    require(f.rows()==guard->faces().rows() && (f.array()==guard->faces().array()).all(),"connectivity_changed");
    auto x=Capture(mesh,geom); checkConstraints(x);
    require(std::isfinite(objective) && objective<=startEnergy+1e-8*std::max(1.,std::abs(startEnergy)),
            "complete_step_objective_increased_or_nonfinite");
    // Catch integration mistakes as well as checking every individual motion.
    require(guard->check(x,x).safe,"accepted_endpoint_not_admissible");
    const double rms=(x-stepStart).norm()/std::sqrt(double(x.rows()));
    require(rms>1e-12*initialEdge,"numerical_stall_not_an_equilibrium_certificate");
    const double a=hard_exclusion::area(x,f),v=hard_exclusion::volume(x,f);
    steps << currentStep << ',' << a << ',' << v << ',' << a/targetArea-1 << ',' << v/targetVolume-1 << ','
          << objective << ',' << objective-startEnergy << ',' << rms << ',' << motionCalls << ',' << limitedCalls << ',' << smallestFraction << std::endl;
    ++acceptedSteps; status("running","last_step_accepted");
}
void FailStep(const MeshPtr& mesh,const GeomPtr& geom,const std::string& reason) {
    if(!Enabled()) return;
    // Retain the geometry at the actual failure before restoring the accepted
    // iteration start. Neither diagnostic is an accepted scientific frame.
    writeDiagnosticMesh("hard_exclusion_failure_state.obj",Capture(mesh,geom),guard->faces());
    if(stepStart.rows()>0) Restore(mesh,geom,stepStart);
    writeDiagnosticMesh("hard_exclusion_restored_state.obj",Capture(mesh,geom),guard->faces());
    status("stopped_inconclusive",reason);
}
void Finish() { if(Enabled()) status("iteration_limit_reached","short_run_complete_not_a_convergence_claim"); }
}}
