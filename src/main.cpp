#include "boundary.hpp"
#include "eigen_ext.hpp"
#include "function.hpp"
#include "integrate_func.hpp"
#include "parameters.hpp"
#include "stiff.hpp"

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/Sparse>

// TYPE DEFINITIONS
typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXi MatI;
typedef Eigen::VectorXi VecI;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;

int main(int args, char *argv[]) {

  // PRELIMINARIES
  int max_iter = 1E+2; // 1E+3;
  double tol_min = 1E-10;
  double tol = 1;

  // EXTRACT DATA
  MatI nodes =
      pear::load_csv<MatI>("/Users/hdeplaen/Documents/KULeuven/Project/"
                           "the_winning_pear/imports/node.csv");
  MatI boundary =
      pear::load_csv<MatI>("/Users/hdeplaen/Documents/KULeuven/Project/"
                           "the_winning_pear/imports/boundary.csv");
  Mat points = pear::load_csv<Mat>("/Users/hdeplaen/Documents/KULeuven/Project/"
                                   "the_winning_pear/imports/points.csv");

  int np = points.rows();
  Vec xp(np);
  Vec yp(np);
  xp = points.col(0);
  yp = points.col(1);

  // INITIALIZE MATRICES
  Mat Ku(np, np);
  Mat Kv(np, np);
  Mat int_F(np, np);
  Vec Bu(np);
  Vec Bv(np);
  Mat Kbu(np, np);
  Mat Kbv(np, np);

  pear::stiff(xp, yp, nodes, Ku, Kv);
  int_F = pear::int_func(xp, yp, nodes);
  pear::boundary_vector(xp, yp, boundary, Bu, Bv, Kbu, Kbv);

  // INITIAL SOLUTION
  Vec Cu(np);
  Vec Cv(np);
  Vec sol(2 * np);
  Vec sol_buff(2 * np);

  Eigen::BiCGSTAB<Mat> solver;

  solver.compute(Ku + (pear::Vmu / pear::Kmu) * int_F - Kbu);
  Cu = solver.solve(Bu);

  /* std::cout << Cu << std::endl; */

  solver.compute(Kv - Kbv);
  Cv = solver.solve(pear::rq * (pear::Vmu / pear::Kmu) * int_F * Cu + Bv);

  std::cout << "Initial solution computed" << std::endl << std::endl;

  // NON-LINEAR SOLUTION
  Vec Ru(np);
  Vec Rv(np);
  Mat RudCu(np, np);
  Mat RudCv(np, np);
  Mat RvdCu(np, np);
  Mat RvdCv(np, np);

  Vec F(2 * np);
  Mat JF(2 * np, 2 * np);

  for (int iter = 0; iter < max_iter; iter++) {
    pear::fun_diff(Cu, Cv, Ru, Rv, RudCu, RudCv, RvdCu, RvdCv);

    F.block(0, 0, np, 1) = Ku * Cu + int_F * Ru - Kbu * Cu - Bu;
    F.block(np, 0, np, 1) = Kv * Cv - int_F * Rv - Kbv * Cv - Bv;

    JF.block(0, 0, np, np) = Ku + int_F * RudCu - Kbu;
    JF.block(0, np, np, np) = int_F * RudCv;
    JF.block(np, 0, np, np) = -int_F * RvdCu;
    JF.block(np, np, np, np) = Kv - int_F * RvdCv - Kbv;

    solver.compute(JF);

    sol_buff.block(0, 0, np, 1) = Cu;
    sol_buff.block(np, 0, np, 1) = Cv;

    sol = sol_buff - solver.solve(F);

    tol = (sol - sol_buff).norm();

    std::cout << "Iteration: " << iter + 1 << std::endl;
    std::cout << "Tolerance: " << tol << std::endl << std::endl;

    if (tol < tol_min) {
      break;
    }

    Cu = sol.block(0, 0, np, 1);
    Cv = sol.block(np, 0, np, 1);
  }

  // EXPORTS
  pear::write_csv<Vec>("/Users/hdeplaen/Documents/KULeuven/Project/"
                       "the_winning_pear/exports/cu.csv",
                       Cu);

  pear::write_csv<Vec>("/Users/hdeplaen/Documents/KULeuven/Project/"
                       "the_winning_pear/exports/cv.csv",
                       Cv);
  return 0;
}
