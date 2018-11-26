#ifndef pear_stiff_hpp
#define pear_stiff_hpp

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace pear {
Eigen::MatrixXd grad_phi(Eigen::Vector3d x, Eigen::Vector3d y);
Eigen::MatrixXd stiff(Eigen::VectorXd &xp, Eigen::VectorXd &yp,
                      Eigen::MatrixXi &node, int boundary_vertices, double Dr,
                      double Dz, double h, double Camb);
} // namespace pear
#endif
