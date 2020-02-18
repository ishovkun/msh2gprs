#ifdef WITH_EIGEN
#include "JacobiPreconditioner.hpp"

namespace discretization {

JacobiPreconditioner::JacobiPreconditioner(const Eigen::SparseMatrix<double,Eigen::RowMajor> & mat)
: _invdiag(mat.cols())
{
  factorize_(mat);
}

void JacobiPreconditioner::factorize_(const Eigen::SparseMatrix<double,Eigen::RowMajor> & mat)
{
 _invdiag.resize(mat.cols());
 for(int j=0; j<mat.outerSize(); ++j)
 {
   Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(mat,j);
   while(it && it.index()!=j) ++it;
   if(it && it.index()==j && it.value()!=double(0))
     _invdiag(j) = 1.0 / it.value();
   else
     _invdiag(j) = 1.0;
 }
}

void JacobiPreconditioner::solve(const Eigen::SparseMatrix<double,Eigen::RowMajor> & A,
                                 const Eigen::VectorXd& b, Eigen::VectorXd& x) const
{
  x = - (2.0/3.0) * _invdiag.array() * (A * b).array() ;
}


}  // end namespace discretization

#endif
