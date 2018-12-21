#include "boundary.hpp"
#include "eigen_ext.hpp"
#include "parameters.hpp"
#include <cassert>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace pear {
void boundary_vector(Vec &xp, Vec &yp, MatI &boundary, Vec &Bu, Vec &Bv,
                     Mat &Kbu, Mat &Kbv) {
  /*
  INPUTS
  xp : vector with points x coordinates
  yp : vector with points y coordinates
  g1 : vector of indices of first points of the gamma1 boundary segment
  g2 : vector of indices of second ponts of the gamma1 boundary segment

  OUTPUTS
  Bu, Bv : boundary vector
  Kbu, Kbv : boundary vector term dependent on Cu and Cv
   */
  int np = xp.rows();
  int nb = boundary.rows();
  int ng = 2 * nb - boundary.col(2).sum();

  VecI g1(ng);
  VecI g2(ng);
  compute_skin(boundary, g1, g2);

  // COMPUTE DISTANCES
  Vec hpx =
      (pear::extract<Vec>(xp, g1) - pear::extract<Vec>(xp, g2)).array().pow(2);
  Vec hpy =
      (pear::extract<Vec>(yp, g1) - pear::extract<Vec>(yp, g2)).array().pow(2);

  Vec h = (hpy + hpx).array().pow(.5);

  // COMPUTE BOUNDARY FUNC
  Vec bu(np);
  Vec bv(np);
  Vec kbu(np);
  Vec kbv(np);
  boundary_func(xp, bu, bv, kbu, kbv);

  // BB
  Vec bbu1 =
      (h.array() *
       (2 * pear::extract<Vec>(bu, g1) + pear::extract<Vec>(bu, g2)).array() /
       6);

  Vec bbu2 =
      (h.array() *
       (pear::extract<Vec>(bu, g1) + 2 * pear::extract<Vec>(bu, g2)).array() /
       6);

  Vec bbv1 =
      (h.array() *
       (2 * pear::extract<Vec>(bv, g1) + pear::extract<Vec>(bv, g2)).array() /
       6);

  Vec bbv2 =
      (h.array() *
       (pear::extract<Vec>(bu, g1) + 2 * pear::extract<Vec>(bu, g2)).array() /
       6);

  // KBB
  Vec kbbu1 =
      (h.array() *
       (3 * pear::extract<Vec>(kbu, g1) + pear::extract<Vec>(kbu, g2)).array() /
       12);

  Vec kbbu2 =
      (h.array() *
       (pear::extract<Vec>(kbu, g1) + 3 * pear::extract<Vec>(kbu, g2)).array() /
       12);

  Vec kbbum =
      (h.array() *
       (pear::extract<Vec>(kbu, g1) + pear::extract<Vec>(kbu, g2)).array() /
       12);

  Vec kbbv1 =
      (h.array() *
       (3 * pear::extract<Vec>(kbv, g1) + pear::extract<Vec>(kbv, g2)).array() /
       12);

  Vec kbbv2 =
      (h.array() *
       (pear::extract<Vec>(kbv, g1) + 3 * pear::extract<Vec>(kbv, g2)).array() /
       12);

  Vec kbbvm =
      (h.array() *
       (pear::extract<Vec>(kbv, g1) + pear::extract<Vec>(kbv, g2)).array() /
       12);

  // FILL BOUNDARY VECTOR
  for (int idxg = 0; idxg < ng; idxg++) {
    Bu(g1(idxg)) += bbu1(idxg);
    Bu(g2(idxg)) += bbu2(idxg);
    Bv(g2(idxg)) += bbv1(idxg);
    Bv(g2(idxg)) += bbv2(idxg);
  }

  // FILL BOUNDARY MATRIX

  for (int idxg = 0; idxg < ng; idxg++) {
    Kbu(g1(idxg), g1(idxg)) += kbbu1(idxg);
    Kbu(g2(idxg), g2(idxg)) += kbbu2(idxg);
    Kbu(g1(idxg), g2(idxg)) += kbbum(idxg);
    Kbu(g2(idxg), g1(idxg)) += kbbum(idxg);

    Kbv(g1(idxg), g1(idxg)) += kbbv1(idxg);
    Kbv(g2(idxg), g2(idxg)) += kbbv2(idxg);
    Kbv(g1(idxg), g2(idxg)) += kbbvm(idxg);
    Kbv(g2(idxg), g1(idxg)) += kbbvm(idxg);
  }
}

void boundary_func(Vec &xp, Vec &bu, Vec &bv, Vec &kbu, Vec &kbv) {
  bu = xp * pear::hu * pear::Cuamb;
  bv = xp * pear::hv * pear::Cvamb;

  kbu = -pear::hu * xp;
  kbv = -pear::hv * xp;
}

void compute_skin(MatI &boundary, VecI &g1, VecI &g2) {
  int nb = boundary.rows();
  int ng = 2 * nb - boundary.col(2).sum();

  int idxg = 0;
  for (int idxb = 0; idxb < nb; idxb++) {
    if (boundary(idxb, 2) == 1) {
      g1(idxg) = boundary(idxb, 0);
      g2(idxg) = boundary(idxb, 1);
      idxg++;
    }
  }
}
} // namespace pear
