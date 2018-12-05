#ifndef pear_boundary_hpp
#define pear_boundary_hpp

#include <Eigen/Core>
#include <Eigen/Sparse>

// TYPE DEFINITIONS
typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXi MatI;
typedef Eigen::VectorXi VecI;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;

namespace pear {
void boundary_vector(Vec &xp, Vec &yp, Vec &g1, Vec &g2, Vec &B, Mat &Kb);
void boundary_func(Vec &xp, Vec &bu, Vec &bv, Vec &kbu, Vec &kbv);
} // namespace pear

#endif
