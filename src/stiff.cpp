#include "eigen_ext.hpp"
#include "parameters.hpp"
#include "stiff.hpp"
#include <cassert>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/SparseCore>

// TYPE DEFINITIONS
typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXi MatI;
typedef Eigen::VectorXi VecI;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;

namespace pear {
void grad_phi(Vec &xp, Vec &yp, Vec &t, Vec &Dphi1, Vec &Dphi2, Vec &Dphi3,
              Vec &T) {
  int np = p.rows();

  VecI t1 = t.col(0);
  VecI t2 = t.col(1);
  VecI t3 = t.col(2);

  Vec r1 = pear::extract<Vec>(xp, t1);
  Vec r2 = pear::extract<Vec>(xp, t2);
  Vec r3 = pear::extract<Vec>(xp, t3);
  Vec z1 = pear::extract<Vec>(yp, t1);
  Vec z2 = pear::extract<Vec>(yp, t2);
  Vec z3 = pear::extract<Vec>(yp, t3);

  Vec r21 = r2 - r1;
  Vec z21 = z2 - z1;
  Vec r32 = r3 - r2;
  Vec z32 = z3 - z2;
  Vec r31 = r3 - r1;
  Vec z31 = z3 - z1;

  // ADDITIONS
  Vec z = (z1 + z2 + z3) / 3;
  Vec r = (r1 + r2 + r3) / 3;

  // AREA OF THE TRIANGLES
  Vec T = (r21 * z31 - z21 * r31) / 2;

  // GRADIENTS OF THE BASIS FUNCTIONS
  Vec Dphi = .5 * [ -2 * z32 + (r32 * z + r2 * z3 - r3 * z2) / r, r32 ];
  Vec Dphi = .5 * [ 2 * z31 - (r31 * z + r1 * z3 - r3 * z1) / r, -r31 ];
  Vec Dphi = .5 * [ -2 * z21 + (r21 * z + r1 * z2 - r2 * z1) / r, r21 ];

  return Dphi;
}

SpMat stiff(Vec &xp, Vec &yp, MatI &node, int boundary_vertices) {
  // PRELIMINARIES
  int np = p.rows();
  int nt = t.rows();

  Eigen::Vector2d Du(np, 2);
  Eigen::Vector2d Dv(np, 2);

  SpMat R(2 * np, 2 * np);

  // GRADIENTS OF BASIS FUNCTIONS
  Mat Dphi = grad_phi(Vec & xp, Vec & yp, Vec & t);

  Du[0] = Dur;
  Du[1] = Duz;
  Dv[0] = Dvr;
  Dv[1] = Dvz;

  // CONSTRUCTION OF THE STIFNESS MATRIX
  for (int idxi = 0; idxi < 3; idxi++) {
    for (int idxj = 0; idxj < idxi; idx2++) {
      for (int idxv = 0; idxv < nt; idxv++) {
        R(t1(idxv, idxi), t2(idxv, idxj)) +=
            Du * (Dphi.row(idxi).array() * Dphi.row(idxj).array()).vector();

        R(t1(idxv, idxi) + np, t2(idxv, idxj) + np) +=
            Dv * (Dphi.row(idxi).array() * Dphi.row(idxj).array()).vector();
        ;
      }
    }
  }

  R += R.transpose();

  return R;
}

} // namespace pear
