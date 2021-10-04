#pragma once
#include "SimData.hpp"
#ifdef WITH_EIGEN
#include <Eigen/Sparse>                     // provides SparseMatrix
#include <Eigen/Dense>                      // provides MatrixXd, VectorXd
#include <vector>

namespace multiscale {

class ShapeFunctionSolver {
 public:
  ShapeFunctionSolver(size_t source, std::vector<size_t> const & bnd_cells, gprs_data::SimData const & data);
  ShapeFunctionSolver(size_t source,
                      std::vector<size_t> const & region,
                      std::vector<size_t> const & bnd,
                      std::vector<double> const & bnd_values,
                      gprs_data::SimData const & data);

  const std::vector<double> & solution() {return _solution;}
  virtual ~ShapeFunctionSolver() = default;

 protected:
  std::tuple<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::VectorXd> build_system_() const;
  void build_laplace_terms_(Eigen::SparseMatrix<double,Eigen::RowMajor>& mat, Eigen::VectorXd & rhs) const;
  void build_special_terms_(Eigen::SparseMatrix<double,Eigen::RowMajor>& mat, Eigen::VectorXd & rhs) const;
  void impose_bc_(Eigen::SparseMatrix<double,Eigen::RowMajor>& mat, Eigen::VectorXd & rhs) const;
  void solve_();

  size_t _source;
  std::vector<size_t> _bnd;
  std::vector<double> _bnd_values;
  gprs_data::SimData const & _data;
  std::vector<double> _solution;
  std::vector<size_t> _region;
  std::vector<size_t> _mapping;
};

}  // end namespace multiscale

#endif
