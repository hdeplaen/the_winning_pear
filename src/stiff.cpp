#include "eigen_ext.hpp"
#include "parameters.hpp"
#include "stiff.hpp"
#include <cassert>
#include <iostream>

using namespace Eigen;

namespace pear {
MatrixXd grad_phi(VectorXd &xp, VectorXd &yp, VectorXd &t1, VectorXd &t2,
                  VectorXd &t3) {
  int np = p.rows();
  MatrixXd Dphi(np, 7);

  VectorXd r1 = pear::extract<VectorXd>(xp, t1);
  VectorXd r2 = pear::extract<VectorXd>(xp, t2);
  VectorXd r3 = pear::extract<VectorXd>(xp, t3);
  VectorXd z1 = pear::extract<VectorXd>(yp, t1);
  VectorXd z2 = pear::extract<VectorXd>(yp, t2);
  VectorXd z3 = pear::extract<VectorXd>(yp, t3);

  VectorXd r21 = r2 - r1;
  VectorXd z21 = z2 - z1;
  VectorXd r32 = r3 - r2;
  VectorXd z32 = z3 - z2;
  VectorXd r31 = r3 - r1;
  VectorXd z31 = z3 - z1;

  // ADDITIONS
  VectorXd z = (z1 + z2 + z3) / 3;
  VectorXd r = (r1 + r2 + r3) / 3;

  // AREA OF THE TRIANGLES
  VectorXd T = (r21 * z31 - z21 * r31) / 2;

  VectorXd dphi1 =
      .5 * [ Dur * (-2 * z32 + (r32 * z + r2 * z3 - r3 * z2) / r), Duz * r32 ];
  VectorXd dphi2 =
      .5 *
      [ Dur * (2 * z31 - (r31 * z + r1 * z3 - r3 * z1) / r), Duz * (-r31) ];
  VectorXd dphi3 =
      .5 *
      [ Dur * (-2 * z21 + (r21 * z + r1 * z2 - r2 * z1) / r), Duz * (r21) ];

  return Dphi;
}

MatrixXd stiff(VectorXd &xp, VectorXd &yp, MatrixXi &node,
               int boundary_vertices) {}
} // namespace pear
