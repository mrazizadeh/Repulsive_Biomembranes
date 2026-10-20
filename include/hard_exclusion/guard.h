#pragma once

#include <Eigen/Core>
#include <string>
#include <vector>

// Independent geometric guard. No energy/force, no changes to connectivity.
// Each accepted displacement is linear; a sequence of separately guarded
// displacements covers descent, Newton correction and barycenter projection.
// PRECONDITION: starting embedding must be independently checked for static
// triangle intersections (CGAL + MeshLab in prepare_cases.py, bound to SHA256).
// Static VF/EE feature queries alone do NOT detect every pre-existing
// edge-through-face penetration. Use the supplied runner and audited inputs.
namespace hard_exclusion {
using Positions = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using Faces = Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>;

struct Options {
    double minimum_separation = 1e-7;
    double tolerance = 1e-9;
    long maximum_ccd_iterations = 100000;
    int maximum_backtracks = 40;
    double normal_dot_floor = 1e-8;
};

struct Report {
    bool safe = false;
    double fraction = 0;
    double possible_toi = 0;
    long face_pairs = 0;
    long vertex_face_queries = 0;
    long edge_edge_queries = 0;
    int backtracks = 0;
    int face_a = -1;
    int face_b = -1;
    std::string reason;
};

class Guard {
public:
    Guard(const Faces& faces, const Options& options);
    Report check(const Positions& start, const Positions& end) const;
    Report filter(const Positions& start, const Positions& proposed) const;
    const Faces& faces() const { return faces_; }
    const Options& options() const { return options_; }
private:
    Faces faces_;
    Options options_;
};

double area(const Positions& x, const Faces& f);
double volume(const Positions& x, const Faces& f);
double mean_edge_length(const Positions& x, const Faces& f);
void validate_closed_manifold(const Positions& x, const Faces& f);
}
