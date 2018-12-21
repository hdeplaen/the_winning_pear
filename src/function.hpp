#ifndef pear_func_hpp
#define pear_func_hpp

#include <Eigen/Core>
#include <Eigen/Sparse>

// TYPE DEFINITIONS
typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXi MatI;
typedef Eigen::VectorXi VecI;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;

namespace pear {
void fun_diff(Vec &Cu, Vec &Cv, Vec &Ru, Vec &Rv, Mat &RudCu, Mat &RudCv,
              Mat &RvdCu, Mat &RvdCv);
} // namespace pear

#endif
