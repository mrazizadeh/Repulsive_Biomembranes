#include "hard_exclusion/linear_solve.h"
#include <iostream>
#include <unsupported/Eigen/IterativeSolvers>
#include <stdexcept>

int main() {
    try {
        Eigen::MatrixXd a=Eigen::MatrixXd::Zero(5,5);
        for(int i=0;i<5;++i) {a(i,i)=i+2; if(i<4) a(i,i+1)=0.7;}
        Eigen::VectorXd b=Eigen::VectorXd::Ones(5),zero=Eigen::VectorXd::Zero(5);
        Eigen::GMRES<Eigen::MatrixXd,Eigen::IdentityPreconditioner> solver;
        solver.setTolerance(1e-4); solver.setMaxIterations(1); solver.compute(a);
        Eigen::VectorXd x=solver.solveWithGuess(b,zero);
        if(!x.allFinite() || solver.info()==Eigen::Success ||
           biorsurfaces::HardExclusion::LinearSolveAcceptable(int(solver.info()),solver.error(),solver.tolerance(),x.allFinite()))
            throw std::runtime_error("Finite but unconverged GMRES result was not rejected");
        std::cout << "PASS finite unconverged GMRES rejected: info=" << int(solver.info()) << " residual=" << solver.error() << '\n';
        solver.setMaxIterations(100); solver.setTolerance(1e-12);
        x=solver.solveWithGuess(b,zero);
        if(!biorsurfaces::HardExclusion::LinearSolveAcceptable(int(solver.info()),solver.error(),solver.tolerance(),x.allFinite()) ||
           (a*x-b).norm()/b.norm()>1e-10)
            throw std::runtime_error("Converged GMRES was rejected or has large independently computed residual");
        std::cout << "PASS converged GMRES accepted; independent residual=" << (a*x-b).norm()/b.norm() << '\n';
        std::cout << "Eigen " << EIGEN_WORLD_VERSION << '.' << EIGEN_MAJOR_VERSION << '.' << EIGEN_MINOR_VERSION << '\n';
    } catch(const std::exception& e) {std::cerr << e.what() << '\n'; return 1;}
    return 0;
}
