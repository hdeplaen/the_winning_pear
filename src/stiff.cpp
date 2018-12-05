#include "eigen_ext.hpp"
#include "parameters.hpp"
#include "stiff.hpp"
#include <cassert>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace pear {
void grad_phi(Vec &xp, Vec &yp, Vec &t, Mat &Dphi1, Mat &Dphi2, Mat &Dphi3,
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
  T = ((r21.array() * z31.array() - z21.array() * r31.array())).vector() / 2;

  // GRADIENTS OF THE BASIS FUNCTIONS
  Dphi1 << -2 * z32.array() + (r32.array() * z + r2.array() * z3.array() -
                               r3.array() * z2.array()) /
                                  r.array(),
      r32.array();

  Dphi2 << 2 * z31.array() -
               (r31.array() * z.array() + r1.array() * z3.array() -
                r3.array() * z1.array()) /
                   r.array(),
      -r31.array();
  Dphi3 << -2 * z21.array() +
               (r21.array() * z.array() + r1.array() * z2.array() -
                r2.array() * z1.array()) /
                   r.array(),
      r21.array();

  Dphi1 = Dphi1 / 2;
  Dphi2 = Dphi2 / 2;
  Dphi3 = Dphi3 / 2;
}

SpMat stiff(Vec &xp, Vec &yp, MatI &node, int boundary_vertices) {
  // PRELIMINARIES
  int np = p.rows();
  int nt = t.rows();

  Eigen::Vector2d Du(np, 2);
  Eigen::Vector2d Dv(np, 2);

  SpMat R(2 * np, 2 * np);

  // GRADIENTS OF BASIS FUNCTIONS
  Mat Dphi1(np, 2);
  Mat Dphi2(np, 2);
  Mat Dphi3(np, 2);
  Vec T(np);

  grad_phi(Vec & xp, Vec & yp, Vec & t, Dphi1, Dphi2, Dphi3, T);

  Du[0] = Dur;
  Du[1] = Duz;
  Dv[0] = Dvr;
  Dv[1] = Dvz;

  // CONSTRUCTION OF THE STIFNESS MATRIX
  for (int idxi = 0; idxi < 3; idxi++) {
    for (int idxj = 0; idxj < idxi; idx2++) {
      for (int idxv = 0; idxv < np; idxv++) {
        R(t1(idxv, idxi), t2(idxv, idxj)) +=
            Du * (Dphi.row(idxi).array() * Dphi.row(idxj).array()).vector();

        R(t1(idxv, idxi) + np, t2(idxv, idxj) + np) +=
            Dv * (Dphi.row(idxi).array() * Dphi.row(idxj).array()).vector();
      }
    }
  }

  R += R.transpose();

  for (int idxi = 0; idxi < 3; idxi++) {
    for (int idxv = 0; idxv < nt; idxv++) {
      R(t1(idxv, idxi), t2(idxv, idxi)) +=
          Du * (Dphi.row(idxi).array() * Dphi.row(idxi).array()).vector();

      R(t1(idxv, idxi) + np, t2(idxv, idxj) + np) +=
          Dv * (Dphi.row(idxi).array() * Dphi.row(idxi).array()).vector();
    }
  }

  return R;
}

} // namespace pear
