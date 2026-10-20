#include "hard_exclusion/guard.h"
#include <Eigen/Geometry>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

using namespace hard_exclusion;
void require(bool condition,const std::string& message) {
    if(!condition) throw std::runtime_error(message);
}
void read_obj(const std::string& path,Positions& x,Faces& f) {
    std::ifstream in(path); if(!in) throw std::runtime_error("Cannot open "+path);
    std::vector<Eigen::Vector3d> vs; std::vector<Eigen::Vector3i> fs;
    std::string line;
    while(std::getline(in,line)) {
        std::istringstream row(line); std::string tag; row>>tag;
        if(tag=="v") {Eigen::Vector3d v; row>>v[0]>>v[1]>>v[2]; require(bool(row),"Bad vertex");vs.push_back(v);}
        if(tag=="f") {
            Eigen::Vector3i face;
            for(int j=0;j<3;++j){std::string token;row>>token;face[j]=std::stoi(token.substr(0,token.find('/')))-1;}
            std::string extra;require(!(row>>extra),"Only triangular OBJ faces are supported");fs.push_back(face);
        }
    }
    x.resize(vs.size(),3);f.resize(fs.size(),3);
    for(size_t i=0;i<vs.size();++i)x.row(i)=vs[i].transpose();
    for(size_t i=0;i<fs.size();++i)f.row(i)=fs[i].transpose();
}
void output(const Report& r,double elapsed) {
    std::cout<<std::setprecision(17)<<"{\"safe\":"<<(r.safe?"true":"false")
        <<",\"fraction\":"<<r.fraction<<",\"reason\":\""<<r.reason
        <<"\",\"face_pairs\":"<<r.face_pairs<<",\"vf_queries\":"<<r.vertex_face_queries
        <<",\"ee_queries\":"<<r.edge_edge_queries<<",\"backtracks\":"<<r.backtracks
        <<",\"seconds\":"<<elapsed<<"}\n";
}
int main(int argc,char** argv) {
    try {
        if(argc>=3) {
            Positions a,b;Faces f,g;read_obj(argv[2],a,f);
            require(a.rows()>0 && f.rows()>0 && a.allFinite() && f.minCoeff()>=0 && f.maxCoeff()<a.rows(),"Invalid OBJ coordinates/indices");
            Options opt;const double h=mean_edge_length(a,f);
            double relative=1e-6;
            if(std::string(argv[1])=="--inspect" && argc==4) relative=std::stod(argv[3]);
            opt.minimum_separation=relative*h;opt.tolerance=1e-8*h;
            Guard guard(f,opt);auto start=std::chrono::steady_clock::now();Report r;
            if(std::string(argv[1])=="--inspect") {validate_closed_manifold(a,f);r=guard.check(a,a);}
            else if(std::string(argv[1])=="--motion" && argc==4) {
                read_obj(argv[3],b,g);require(g.rows()==f.rows() && (g.array()==f.array()).all(),"Connectivity differs");
                r=guard.check(a,b);
            } else throw std::runtime_error("Usage: test_guard [--inspect mesh.obj [relative_clearance] | --motion start.obj end.obj]");
            output(r,std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());return r.safe?0:3;
        }
        Options o;o.minimum_separation=1e-7;o.tolerance=1e-9;
        Faces f(2,3);f<<0,1,2,3,4,5;
        Positions a(6,3);a<<-1,-1,0, 1,-1,0, 0,1,0, -.1,-.1,1, .1,-.1,1, 0,.1,1;
        Guard guard(f,o);require(guard.check(a,a).safe,"Safe initial triangles rejected");
        Positions b=a;b.bottomRows(3).col(2).array()-=2;
        require(guard.check(b,b).safe,"Tunneling endpoint should be clear");
        require(!guard.check(a,b).safe,"Missed tunneling through face");
        Report limited=guard.filter(a,b);require(limited.safe && limited.fraction>0 && limited.fraction<.5,"Failed to shorten crossing motion");
        require(guard.check(a,Positions(a+limited.fraction*(b-a))).safe,"Shortened step is unsafe");
        std::cout<<"PASS vertex-face tunneling and rechecked step reduction\n";
        b=a;b.rowwise()+=Eigen::RowVector3d(.2,.3,.4);
        require(guard.check(a,b).safe,"Safe rigid translation rejected");
        b=a;b.bottomRows(3).col(2).array()+=1;
        require(guard.filter(a,b).fraction==1,"Safe separation was modified");
        std::cout<<"PASS safe rigid motion and separating motion\n";
        Positions ee(6,3);ee<<-1,0,0, 1,0,0, -1,-.1,0, 0,-1,1, 0,1,1, 1,2,1;
        b=ee;b.bottomRows(3).col(2).array()-=2;
        Report er=guard.check(ee,b);
        require(!er.safe && er.reason=="edge_edge_contact_or_uncertainty","Missed edge-edge crossing without vertex-face event");
        std::cout<<"PASS edge-edge interior crossing\n";
        b=a;b.bottomRows(3).col(2).setZero();
        require(!guard.check(b,b).safe,"Intersecting input passed");
        std::cout<<"PASS initially coplanar overlap rejected (general static audit is external)\n";
        Faces adjacent(2,3);adjacent<<0,1,2,1,0,3;
        Positions av(4,3);av<<0,0,0,1,0,0,0,1,0,0,-1,0;
        Guard ag(adjacent,o);require(ag.check(av,av).safe,"Legitimate shared edge rejected");
        Positions fold=av;fold(3,1)=1;
        require(!ag.check(av,fold).safe,"Adjacent fold through shared edge missed");
        std::cout<<"PASS legitimate adjacency and invalid adjacent fold\n";
        Faces sf(1,3);sf<<0,1,2;Positions tri=av.topRows(3),flip=tri;
        flip.row(1)*=-1;flip.row(2)*=-1;
        require(!Guard(sf,o).check(tri,flip).safe,"Intermediate triangle collapse missed");
        std::cout<<"PASS intermediate triangle collapse with nondegenerate endpoints\n";
        // A safe descent followed by an unsafe constraint correction must be
        // checked as two separate motions; checking descent alone is insufficient.
        Positions descent=a;descent.bottomRows(3).col(2).array()-=.1;
        require(guard.check(a,descent).safe,"Safe descent rejected");
        Positions corrected=descent;corrected.bottomRows(3).col(2).array()-=2;
        require(!guard.check(descent,corrected).safe,"Constraint-only crossing missed");
        Report cr=guard.filter(descent,corrected);
        require(cr.safe && cr.fraction<.5,"Constraint correction not shortened");
        std::cout<<"PASS crossing introduced only by constraint correction\n";
        Faces many(30,3); Positions mv(90,3);
        for(int i=0;i<30;++i) {
            many.row(i)<<3*i,3*i+1,3*i+2;
            mv.middleRows(3*i,3)=a.topRows(3);
            mv.middleRows(3*i,3).col(0).array()+=10*(i%15);
            mv.middleRows(3*i,3).col(2).array()+=(i/15)*2;
        }
        Guard big(many,o); require(big.check(mv,mv).safe,"BVH safe separated faces rejected");
        Positions moved=mv; moved.bottomRows(3).col(2).array()-=4;
        require(big.check(moved,moved).safe,"BVH tunneling endpoint not clear");
        require(!big.check(mv,moved).safe,"BVH missed isolated last-face crossing");
        std::cout<<"PASS BVH traversal with isolated crossing among 30 faces\n";
        Positions pinched(7,3);pinched<<0,0,0,1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1;
        Faces shells(8,3);shells<<0,2,1,0,1,3,0,3,2,1,2,3,0,4,5,0,6,4,0,5,6,4,6,5;
        validate_closed_manifold(Positions(pinched.topRows(4)),Faces(shells.topRows(4)));
        bool pinchedRejected=false;
        try {validate_closed_manifold(pinched,shells);} catch(const std::exception&) {pinchedRejected=true;}
        require(pinchedRejected,"Pinched vertex with two closed fans passed topology check");
        std::cout<<"PASS valid closed tetrahedron and pinched-vertex rejection\n";
        b=a;b(0,0)=std::numeric_limits<double>::quiet_NaN();
        require(!guard.check(a,b).safe,"Nonfinite input accepted");
        Options coarse=o;coarse.minimum_separation=1e-6;
        b=a;b.bottomRows(3).col(2).array()-=2;
        require(!Guard(f,coarse).check(a,b).safe,"Crossing missed at second clearance");
        std::cout<<"PASS nonfinite rejection and second separation tolerance\nALL GUARD TESTS PASSED\n";
    }catch(const std::exception& e){std::cerr<<"TEST FAILED: "<<e.what()<<"\n";return 1;}
    return 0;
}
