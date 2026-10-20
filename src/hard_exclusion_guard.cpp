#include "hard_exclusion/guard.h"
#include <tight_inclusion/ccd.hpp>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

namespace hard_exclusion {
namespace {
using V = Eigen::Vector3d;
V point(const Positions& x, int i) { return x.row(i).transpose(); }

struct Box {
    V lo = V::Constant(std::numeric_limits<double>::infinity());
    V hi = V::Constant(-std::numeric_limits<double>::infinity());
    void add(const V& x) { lo = lo.cwiseMin(x); hi = hi.cwiseMax(x); }
    void add(const Box& b) { add(b.lo); add(b.hi); }
    bool overlaps(const Box& b, double pad) const {
        return !((lo.array() > b.hi.array() + pad).any() ||
                 (b.lo.array() > hi.array() + pad).any());
    }
};
struct Node { Box box; int begin, end, left = -1, right = -1; };

class Query {
public:
    const Positions& x;
    const Positions& y;
    const Faces& f;
    const Options& opt;
    std::vector<Box> boxes;
    std::vector<int> order;
    std::vector<Node> nodes;
    double pad;
    Report result;

    Query(const Positions& a, const Positions& b, const Faces& faces, const Options& options)
        : x(a), y(b), f(faces), opt(options) {
        const double extent = std::max(1., std::max(a.cwiseAbs().maxCoeff(), b.cwiseAbs().maxCoeff()));
        pad = 2 * (opt.minimum_separation + opt.tolerance) +
              128 * std::numeric_limits<double>::epsilon() * extent;
        order.resize(f.rows()); std::iota(order.begin(), order.end(), 0);
        boxes.resize(f.rows()); nodes.reserve(2 * f.rows());
        for (int i = 0; i < f.rows(); ++i)
            for (int j = 0; j < 3; ++j) {
                boxes[i].add(point(x, f(i,j))); boxes[i].add(point(y, f(i,j)));
            }
    }
    int build(int begin, int end) {
        const int id = static_cast<int>(nodes.size()); nodes.emplace_back();
        nodes[id].begin = begin; nodes[id].end = end;
        for (int k = begin; k < end; ++k) nodes[id].box.add(boxes[order[k]]);
        if (end - begin > 4) {
            Eigen::Index axis; (nodes[id].box.hi - nodes[id].box.lo).maxCoeff(&axis);
            const int mid = (begin + end) / 2;
            std::nth_element(order.begin()+begin, order.begin()+mid, order.begin()+end,
                [&](int a, int b) { return boxes[a].lo[axis]+boxes[a].hi[axis] < boxes[b].lo[axis]+boxes[b].hi[axis]; });
            const int l = build(begin, mid), r = build(mid, end);
            nodes[id].left = l; nodes[id].right = r;
        }
        return id;
    }
    Box primitive_box(const int* ids, int n) const {
        Box b;
        for (int i = 0; i < n; ++i) { b.add(point(x,ids[i])); b.add(point(y,ids[i])); }
        return b;
    }
    bool vertex_face(int v, const int* t) {
        if (v == t[0] || v == t[1] || v == t[2]) return false;
        if (!primitive_box(&v,1).overlaps(primitive_box(t,3),pad)) return false;
        ++result.vertex_face_queries;
        double toi = 0, achieved = 0;
        const bool hit = ticcd::vertexFaceCCD(point(x,v),point(x,t[0]),point(x,t[1]),point(x,t[2]),
            point(y,v),point(y,t[0]),point(y,t[1]),point(y,t[2]),Eigen::Array3d::Constant(-1),
            opt.minimum_separation,toi,opt.tolerance,1.,opt.maximum_ccd_iterations,achieved,false,
            ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);
        if (hit) { result.possible_toi = toi; result.reason = "vertex_face_contact_or_uncertainty"; }
        return hit;
    }
    bool edge_edge(const int* a, const int* b) {
        if (a[0]==b[0] || a[0]==b[1] || a[1]==b[0] || a[1]==b[1]) return false;
        if (!primitive_box(a,2).overlaps(primitive_box(b,2),pad)) return false;
        ++result.edge_edge_queries;
        double toi = 0, achieved = 0;
        const bool hit = ticcd::edgeEdgeCCD(point(x,a[0]),point(x,a[1]),point(x,b[0]),point(x,b[1]),
            point(y,a[0]),point(y,a[1]),point(y,b[0]),point(y,b[1]),Eigen::Array3d::Constant(-1),
            opt.minimum_separation,toi,opt.tolerance,1.,opt.maximum_ccd_iterations,achieved,false,
            ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);
        if (hit) { result.possible_toi = toi; result.reason = "edge_edge_contact_or_uncertainty"; }
        return hit;
    }
    bool triangles(int a, int b) {
        if (!boxes[a].overlaps(boxes[b],pad)) return false;
        ++result.face_pairs;
        int fa[3] = {f(a,0),f(a,1),f(a,2)}, fb[3] = {f(b,0),f(b,1),f(b,2)};
        // Do not discard a whole triangle pair just because it shares a vertex.
        // Only intrinsically incident VF/EE primitives are omitted.
        for (int i = 0; i < 3; ++i) {
            if (vertex_face(fa[i],fb) || vertex_face(fb[i],fa)) {
                result.face_a=a; result.face_b=b; return true;
            }
        }
        for (int i=0;i<3;++i) for (int j=0;j<3;++j) {
            int ea[2]={fa[i],fa[(i+1)%3]}, eb[2]={fb[j],fb[(j+1)%3]};
            if (edge_edge(ea,eb)) { result.face_a=a; result.face_b=b; return true; }
        }
        return false;
    }
    bool visit(int a, int b) {
        const Node& na=nodes[a]; const Node& nb=nodes[b];
        if (!na.box.overlaps(nb.box,pad)) return false;
        if (a==b) {
            if (na.left<0) {
                for (int i=na.begin;i<na.end;++i) for (int j=i+1;j<na.end;++j)
                    if (triangles(order[i],order[j])) return true;
                return false;
            }
            return visit(na.left,na.left) || visit(na.left,na.right) || visit(na.right,na.right);
        }
        if (na.left<0 && nb.left<0) {
            for (int i=na.begin;i<na.end;++i) for (int j=nb.begin;j<nb.end;++j)
                if (triangles(order[i],order[j])) return true;
            return false;
        }
        if (nb.left<0 || (na.left>=0 && na.end-na.begin >= nb.end-nb.begin))
            return visit(na.left,b) || visit(na.right,b);
        return visit(a,nb.left) || visit(a,nb.right);
    }
};
}

Guard::Guard(const Faces& f, const Options& o) : faces_(f), options_(o) {
    if (f.rows()==0 || !std::isfinite(o.minimum_separation) || o.minimum_separation<=0 || o.minimum_separation>=1 ||
        !std::isfinite(o.tolerance) || o.tolerance<=0 || o.tolerance>=o.minimum_separation ||
        o.maximum_ccd_iterations<=0 || o.maximum_backtracks<1 || !std::isfinite(o.normal_dot_floor) ||
        o.normal_dot_floor<=0 || o.normal_dot_floor>=1)
        throw std::runtime_error("Invalid hard-exclusion options");
}

Report Guard::check(const Positions& a, const Positions& b) const {
    Report r;
    if (a.rows()==0 || a.rows()!=b.rows() || !a.allFinite() || !b.allFinite() ||
        faces_.minCoeff()<0 || faces_.maxCoeff()>=a.rows()) {
        r.reason="invalid_coordinates_or_indices"; return r;
    }
    // A moving face normal is quadratic in t. Its dot with its initial normal
    // has an analytically computable minimum. Positive dot throughout prevents
    // degeneration/inversion, including a collapse hidden between endpoints.
    for (int i=0;i<faces_.rows();++i) {
        int i0=faces_(i,0),i1=faces_(i,1),i2=faces_(i,2);
        V u=point(a,i1)-point(a,i0), v=point(a,i2)-point(a,i0);
        V du=point(b,i1)-point(b,i0)-u, dv=point(b,i2)-point(b,i0)-v;
        V n=u.cross(v);
        const double q0=n.squaredNorm(), q1=n.dot(du.cross(v)+u.cross(dv)), q2=n.dot(du.cross(dv));
        double minimum=std::min(q0,q0+q1+q2);
        if (q2>0) { double t=-q1/(2*q2); if (t>0 && t<1) minimum=std::min(minimum,q0+t*(q1+t*q2)); }
        if (!std::isfinite(minimum) || q0<=0 || minimum <= options_.normal_dot_floor*q0) {
            r.reason="triangle_degeneracy_or_large_normal_rotation"; r.face_a=i; return r;
        }
    }
    Query q(a,b,faces_,options_); const int root=q.build(0,faces_.rows());
    if (q.visit(root,root)) return q.result;
    q.result.safe=true; q.result.fraction=1; q.result.possible_toi=std::numeric_limits<double>::infinity();
    q.result.reason="certified_motion"; return q.result;
}

Report Guard::filter(const Positions& a, const Positions& proposed) const {
    double fraction=1; Report total;
    for (int k=0;k<options_.maximum_backtracks;++k) {
        Positions candidate=a+fraction*(proposed-a);
        Report r=check(a,candidate);
        total.face_pairs+=r.face_pairs; total.vertex_face_queries+=r.vertex_face_queries;
        total.edge_edge_queries+=r.edge_edge_queries; total.backtracks=k;
        total.reason=r.reason; total.face_a=r.face_a; total.face_b=r.face_b;
        total.possible_toi=r.possible_toi;
        if (r.safe) { total.safe=true; total.fraction=fraction; return total; }
        // Always rerun CCD on the shortened motion. Do not assume that a
        // returned approximate TOI by itself certifies the accepted interval.
        double reduction=.5;
        if (std::isfinite(r.possible_toi) && r.possible_toi>0)
            reduction=std::min(.5,.8*r.possible_toi);
        fraction*=reduction;
        if (fraction<1e-12) break;
    }
    total.safe=false; total.fraction=0; total.reason="no_admissible_motion:"+total.reason; return total;
}

double area(const Positions& x,const Faces& f) {
    double a=0; for(int i=0;i<f.rows();++i)
        a+=(point(x,f(i,1))-point(x,f(i,0))).cross(point(x,f(i,2))-point(x,f(i,0))).norm()/2;
    return a;
}
double volume(const Positions& x,const Faces& f) {
    double v=0; for(int i=0;i<f.rows();++i)
        v+=point(x,f(i,0)).dot(point(x,f(i,1)).cross(point(x,f(i,2))))/6;
    return v;
}
double mean_edge_length(const Positions& x,const Faces& f) {
    double length=0; for(int i=0;i<f.rows();++i) for(int j=0;j<3;++j)
        length+=(point(x,f(i,j))-point(x,f(i,(j+1)%3))).norm();
    return length/(3*f.rows());
}
void validate_closed_manifold(const Positions& x,const Faces& f) {
    if (f.rows()==0 || x.rows()==0 || !x.allFinite() || f.minCoeff()<0 || f.maxCoeff()>=x.rows())
        throw std::runtime_error("Invalid input mesh indices/coordinates");
    std::map<std::pair<int,int>,std::pair<int,int>> edges;
    std::vector<std::vector<int>> neighbors(x.rows());
    std::vector<std::map<int,std::vector<int>>> links(x.rows());
    for(int i=0;i<f.rows();++i) for(int j=0;j<3;++j) {
        int a=f(i,j),b=f(i,(j+1)%3);
        if(a==b) throw std::runtime_error("Degenerate input edge");
        auto& entry=edges[std::minmax(a,b)]; ++entry.first; entry.second+=(a<b?1:-1);
        neighbors[a].push_back(b); neighbors[b].push_back(a);
        int c=f(i,(j+2)%3);
        links[a][b].push_back(c); links[a][c].push_back(b);
    }
    for(const auto& e:edges) if(e.second.first!=2 || e.second.second!=0)
        throw std::runtime_error("Input must be a closed consistently oriented two-manifold");
    // Edge incidence alone misses two closed fans pinched at one vertex.
    // Every vertex link must be one connected cycle.
    for(const auto& link:links) {
        if(link.empty()) throw std::runtime_error("Unused input vertex");
        for(const auto& entry:link) if(entry.second.size()!=2)
            throw std::runtime_error("Input vertex link is not a cycle");
        std::vector<int> walk(1,link.begin()->first); std::set<int> visited(walk.begin(),walk.end());
        for(size_t k=0;k<walk.size();++k) for(int n:link.at(walk[k]))
            if(visited.insert(n).second) walk.push_back(n);
        if(visited.size()!=link.size()) throw std::runtime_error("Pinched nonmanifold input vertex");
    }
    std::vector<bool> seen(x.rows(),false); std::vector<int> todo(1,0); seen[0]=true;
    for(size_t i=0;i<todo.size();++i) for(int n:neighbors[todo[i]]) if(!seen[n]) {seen[n]=true;todo.push_back(n);}
    if(todo.size()!=static_cast<size_t>(x.rows()) || volume(x,f)<=0)
        throw std::runtime_error("Input must be connected, have no unused vertices and have positive oriented volume");
}
}
