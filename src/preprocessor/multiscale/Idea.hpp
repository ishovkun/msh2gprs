#pragma once
#include "mesh/Mesh.hpp"
#include "SimData.hpp"
#include <Eigen/Sparse>                     // provides SparseMatrix
#include <Eigen/Dense>                      // provides MatrixXd, VectorXd
#include <vector>

namespace multiscale {

class Idea {
 public:
  Idea(mesh::Mesh const & grid, gprs_data::SimData & data);
  virtual ~Idea() = default;

 protected:
  size_t find_center_cell_() const;
  std::vector<size_t> find_boundary_cells_() const;
  std::tuple<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::VectorXd> build_system_() const;
  void build_laplace_terms_(Eigen::SparseMatrix<double,Eigen::RowMajor>& mat, Eigen::VectorXd & rhs) const;
  void build_special_terms_(Eigen::SparseMatrix<double,Eigen::RowMajor>& mat, Eigen::VectorXd & rhs) const;
  void impose_bc_(Eigen::SparseMatrix<double,Eigen::RowMajor>& mat, Eigen::VectorXd & rhs) const;

  mesh::Mesh const & _grid;
  gprs_data::SimData & _data;
  size_t _source;
  std::vector<size_t> _boundary_cells;
};

}  // end namespace multiscale
