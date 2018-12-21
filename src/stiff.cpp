#include "eigen_ext.hpp"
#include "parameters.hpp"
#include "stiff.hpp"
#include <cassert>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace pear {
void grad_phi(Vec &xp, Vec &yp, double &T, Vec &Dphi2, Vec &Dphi3) {

  Dphi2 << yp(1) - yp(2), yp(2) - yp(0), yp(0) - yp(1);
  Dphi3 << xp(2) - xp(1), xp(0) - xp(2), xp(1) - xp(0);

  T = xp(1) * yp(2) + xp(0) * yp(1) + xp(2) * yp(0) - xp(1) * yp(0) -
      xp(0) * yp(2) - xp(2) * yp(1);
}

Mat stiff_block(Vec &xp, Vec &yp, double Dr, double Dz) {

  Mat K_block(3, 3);

  Vec Dphi2(3);
  Vec Dphi3(3);
  double T = 1;

  grad_phi(xp, yp, T, Dphi2, Dphi3);

  for (int idx1 = 0; idx1 < 3; idx1++) {
    for (int idx2 = 0; idx2 < 3; idx2++) {
      K_block(idx1, idx2) =
          (xp(0) + xp(1) + xp(2)) *
          (Dr * Dphi2(idx1) * Dphi2(idx2) + Dz * Dphi3(idx1) * Dphi3(idx2)) /
          12 / T;
    }
  }
  return K_block;
}

void stiff(Vec &xp, Vec &yp, MatI &t, Mat &Ku, Mat &Kv) {
  // PRELIMINARIES
  int np = xp.rows();
  int nt = t.rows();

  Mat Ku_block(3, 3);
  Mat Kv_block(3, 3);

  VecI t_loc(3);
  Vec xp_loc(3);
  Vec yp_loc(3);

  for (int idxm = 0; idxm < nt; idxm++) {
    t_loc = t.row(idxm);
    xp_loc = pear::extract<Vec>(xp, t_loc);
    yp_loc = pear::extract<Vec>(yp, t_loc);

    Ku_block = stiff_block(xp_loc, yp_loc, pear::Dur, pear::Duz);
    Kv_block = stiff_block(xp_loc, yp_loc, pear::Dvr, pear::Dvz);
    for (int idx1 = 0; idx1 < 3; idx1++) {
      for (int idx2 = 0; idx2 < 3; idx2++) {
        Ku(t_loc(idx1), t_loc(idx2)) += Ku_block(idx1, idx2);
        Kv(t_loc(idx1), t_loc(idx2)) += Kv_block(idx1, idx2);
      }
    }
  }
}

} // namespace pear
