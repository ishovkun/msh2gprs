#pragma once
#include <Eigen/Sparse>  // provides SparseMatrix
#include <Eigen/Dense>  // provides MatrixXd, VectorXd

namespace discretization {

class JacobiPreconditioner
{
 public:
  JacobiPreconditioner(const Eigen::SparseMatrix<double,Eigen::RowMajor> & mat);
  void solve(const Eigen::SparseMatrix<double,Eigen::RowMajor> & mat,
             const Eigen::VectorXd& b, Eigen::VectorXd& x) const;

 protected:
  void factorize_(const Eigen::SparseMatrix<double,Eigen::RowMajor> & mat);

 private:
  Eigen::VectorXd _invdiag;
};


}  // end namespace discretization
