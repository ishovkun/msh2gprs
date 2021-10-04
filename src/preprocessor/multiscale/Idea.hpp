#pragma once
#include "mesh/Mesh.hpp"
#include "SimData.hpp"
#ifdef WITH_EIGEN
#include <Eigen/Sparse>                     // provides SparseMatrix
#include <Eigen/Dense>                      // provides MatrixXd, VectorXd
#endif
#include <vector>

namespace multiscale {

class Idea {
 public:
  Idea(mesh::Mesh const & grid, gprs_data::SimData & data);
  virtual ~Idea() = default;

 protected:
  size_t find_center_cell_() const;
  std::vector<size_t> find_boundary_cells_() const;
  std::vector<size_t> find_boundary_cells_(std::vector<int> const & bnd_markers );
  void debug_save_solution_(std::string const & fname, std::vector<double>const&) const;

  mesh::Mesh const & _grid;
  gprs_data::SimData & _data;
};

}  // end namespace multiscale
