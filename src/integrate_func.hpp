#ifndef pear_int_hpp
#define pear_int_hpp

#include <Eigen/Core>
#include <Eigen/Sparse>

// TYPE DEFINITIONS
typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXi MatI;
typedef Eigen::VectorXi VecI;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;

namespace pear {
Mat int_func(Vec &xp, Vec &yp, MatI &t);
Mat int_func_block(Vec &xp, Vec &yp);
double triangle_area(Vec &xp, Vec &yp);
} // namespace pear

#endif
