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
void grad_phi(Vec &xp, Vec &yp, double &T, Vec &Dphi2, Vec &Dphi3);
Mat stiff_block(Vec &xp, Vec &yp, double Dr, double Dz);
void stiff(Vec &xp, Vec &yp, MatI &t, Mat &Ku, Mat &Kv);
} // namespace pear

#endif
