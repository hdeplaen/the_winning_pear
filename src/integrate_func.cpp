#include "eigen_ext.hpp"
#include "integrate_func.hpp"
#include "parameters.hpp"
#include <cassert>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace pear {
double triangle_area(Vec &xp, Vec &yp) {
  double area;
  area = xp(1) * yp(2) + xp(0) * yp(1) + xp(2) * yp(0) - xp(1) * yp(0) -
         xp(0) * yp(2) - xp(2) * yp(1);
  return area;
}

Mat int_func_block(Vec &xp, Vec &yp) {
  double area = triangle_area(xp, yp);
  Mat int_func_block(3, 3);

  double b11 = 6 * xp(0) + 2 * xp(1) + 2 * xp(2);
  double b12 = 2 * xp(0) + 2 * xp(1) + xp(2);
  double b13 = 2 * xp(0) + xp(1) + 2 * xp(2);

  double b21 = 2 * xp(0) + 2 * xp(1) + xp(2);
  double b22 = 2 * xp(0) + 6 * xp(1) + 2 * xp(2);
  double b23 = xp(0) + 2 * xp(1) + 2 * xp(2);

  double b31 = 2 * xp(0) + xp(1) + 2 * xp(2);
  double b32 = xp(0) + 2 * xp(1) + 2 * xp(2);
  double b33 = 2 * xp(0) + 2 * xp(1) + 6 * xp(2);

  int_func_block << b11, b12, b13, b21, b22, b23, b31, b32, b33;
  int_func_block = int_func_block * area / 60;

  return int_func_block;
}

Mat int_func(Vec &xp, Vec &yp, MatI &t) {

  int np = xp.rows();
  int nt = t.rows();

  Mat int_F(np, np);
  Mat int_F_block(3, 3);

  VecI t_loc(3);
  Vec xp_loc(3);
  Vec yp_loc(3);

  for (int idxm = 0; idxm < nt; idxm++) {
    t_loc = t.row(idxm);
    xp_loc = pear::extract<Vec>(xp, t_loc);
    yp_loc = pear::extract<Vec>(yp, t_loc);

    int_F_block = int_func_block(xp_loc, yp_loc);
    for (int idx1 = 0; idx1 < 3; idx1++) {
      for (int idx2 = 0; idx2 < 3; idx2++) {
        int_F(t_loc(idx1), t_loc(idx2)) += int_F_block(idx1, idx2);
      }
    }
  }
  return int_F;
}

} // namespace pear
