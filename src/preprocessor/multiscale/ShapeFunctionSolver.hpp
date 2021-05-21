#pragma once
#include "mesh/Mesh.hpp"
#include "SimData.hpp"
#include <Eigen/Sparse>                     // provides SparseMatrix
#include <Eigen/Dense>                      // provides MatrixXd, VectorXd
#include <vector>

namespace multiscale {

class ShapeFunctionSolver {
 public:
  ShapeFunctionSolver(size_t source, std::vector<size_t> const & bnd_cells,
                      mesh::Mesh const & grid, gprs_data::SimData & data);
  const std::vector<double> & solution() {return _solution;}
  virtual ~ShapeFunctionSolver() = default;

 protected:
  std::tuple<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::VectorXd> build_system_() const;
  void build_laplace_terms_(Eigen::SparseMatrix<double,Eigen::RowMajor>& mat, Eigen::VectorXd & rhs) const;
  void build_special_terms_(Eigen::SparseMatrix<double,Eigen::RowMajor>& mat, Eigen::VectorXd & rhs) const;
  void impose_bc_(Eigen::SparseMatrix<double,Eigen::RowMajor>& mat, Eigen::VectorXd & rhs) const;

  size_t _source;
  std::vector<size_t> _boundary_cells;
  mesh::Mesh const & _grid;
  gprs_data::SimData & _data;
  std::vector<double> _solution;
};

}  // end namespace multiscale
