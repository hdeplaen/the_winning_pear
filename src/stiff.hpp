#ifndef pear_stiff_hpp
#define pear_stiff_hpp

#include <Eigen/Core>
#include <Eigen/Sparse>

// TYPE DEFINITIONS
typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXi MatI;
typedef Eigen::VectorXi VecI;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;

namespace pear {
void grad_phi(Vec &xp, Vec &yp, Vec &t, Mat &Dphi1, Mat &Dphi2, Mat &Dphi3,
              Vec &T);
SpMat stiff(Vec &xp, Vec &yp, MatI &node, int boundary_vertices);
} // namespace pear

#endif
