#include "eigen_ext.hpp"
#include "function.hpp"
#include "parameters.hpp"
#include <cassert>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace pear {
void fun_diff(Vec &Cu, Vec &Cv, Vec &Ru, Vec &Rv, Mat &RudCu, Mat &RudCv,
              Mat &RvdCu, Mat &RvdCv) {

  int np = Cu.rows();

  // FUNCTIONS
  Ru = pear::Vmu * Cu.array() /
       ((pear::Kmu + Cu.array()) * (1 + Cv.array() / pear::Kmv));

  Rv = pear::rq * Ru.array() + pear::Vmfv / (1 + Cu.array() / pear::Kmfu);

  // DERIVATIVES
  Vec RudCu_diag(np);
  Vec RudCVdiag(np);
  Vec RvdCu_diag(np);
  Vec RvdCVdiag(np);

  RudCu_diag = (pear::Kmu * pear::Kmv * pear::Vmu) /
               ((Cu.array() + pear::Kmu).pow(2) * (Cv.array() + pear::Kmv));

  RudCVdiag = -(Cu.array() * pear::Kmv * pear::Vmu) /
              ((Cu.array() + pear::Kmu) * (Cv.array() + pear::Kmv).pow(2));

  RvdCu_diag = (pear::Kmv * pear::Vmu * pear::rq) /
                   ((Cu.array() + pear::Kmu) * (Cv.array() + pear::Kmv)) -
               (pear::Kmfu * pear::Vmfv) / (Cu.array() + pear::Kmfu).pow(2) -
               (Cu.array() * pear::Kmv * pear::Vmu * pear::rq) /
                   ((Cu.array() + pear::Kmu).pow(2) * (Cv.array() + pear::Kmv));

  RvdCVdiag = -(Cu.array() * pear::Kmv * pear::Vmu * pear::rq) /
              ((Cu.array() + pear::Kmu) * (Cv.array() + pear::Kmv).pow(2));

  RudCu = RudCu_diag.asDiagonal();
  RvdCu = RudCu_diag.asDiagonal();
  RudCv = RudCu_diag.asDiagonal();
  RvdCv = RudCu_diag.asDiagonal();
}

} // namespace pear
