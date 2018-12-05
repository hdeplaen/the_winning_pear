#include "eigen_ext.hpp"
#include "parameters.hpp"
#include "stiff.hpp"
#include <cassert>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace pear {
void boundary_vector(Vec &xp, Vec &yp, Vec &g1, Vec &g2, Vec &B, Mat &Kb) {
  /*
  INPUTS
  xp : vector with points x coordinates
  yp : vector with points y coordinates
  g1 : vector of indices of first points of the gamma1 boundary segment
  g2 : vector of indices of second ponts of the gamma1 boundary segment

  OUTPUTS
  B : boundary vector
  Kb : boundary vector term dependent on Cu and Cv
   */
  int np = p.rows();
  int ng = g1.rows();

  // COMPUTE DISTANCES
  Vec hpx = (pear::extract<Vec>(xp, g1) - pear::extract<Vec>(xp, g2))
                .array()
                .pow(2)
                .vector();
  Vec hpy = (pear::extract<Vec>(yp, g1) - pear::extract<Vec>(yp, g2))
                .array()
                .pow(2)
                .vector();

  Vec h = (hpy + hpx).array().pow(.5).vector();

  // COMPUTE BOUNDARY FUNC
  Vec bu(np);
  Vec bv(np);
  Vec kbu(np);
  Vec kbv(np);
  boundary_func(xp, bu, bv, kbu, kbv);

  // VALUE OF SEGMENTS
  Vec bbu =
      (h.array() *
       (pear::extract<Vec>(bu, g1) - pear::extract<Vec>(bu, g2)).array() / 2)
          .vector();

  Vec bbv =
      (h.array() *
       (pear::extract<Vec>(bv, g1) - pear::extract<Vec>(bv, g2)).array() / 2)
          .vector();

  Vec kbbu =
      (h.array() *
       (pear::extract<Vec>(kbu, g1) - pear::extract<Vec>(kbu, g2)).array() / 2)
          .vector();

  Vec kbbv =
      (h.array() *
       (pear::extract<Vec>(kbv, g1) - pear::extract<Vec>(kbv, g2)).array() / 2)
          .vector();

  // FILL BOUNDARY VECTOR
  for (int idxg = 0; idxv < ng; idxv++) {
    B(g1(idxg)) += bbu(idxg);
    B(g2(idxg)) += bbu(idxg);
    B(g2(idxg) + np) += bbv(idxg);
    B(g2(idxg) + np) += bbv(idxg);
  }

  // FILL BOUNDARY MATRIX
  for (int idxg = 0; idxv < ng; idxv++) {
    Kb(g1(idxg), g1(idxg)) += kbbu(idxg);
    Kb(g2(idxg), g2(idxg)) += kbbu(idxg);
    Kb(g1(idxg) + np, g1(idxg) + np) += kbbv(idxg);
    Kb(g2(idxg) + np, g2(idxg) + np) += kbbv(idxg);
  }
}

void boundary_func(Vec &xp, Vec &bu, Vec &bv, Vec &kbu, Vec &kbv) {
  bu = xp * pear::hu * pear::Cuamb;
  bv = xp * pear::hv * pear::Cvamb;

  kbu = -pear::hu * xp;
  kbv = -pear::hv * xp;
}

} // namespace pear
