#pragma once

#include "mesh/Cell.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include <Eigen/Sparse>

namespace discretization {

class DFEMElement
{
 public:
  DFEMElement(const mesh::Cell & cell);
 protected:
  void build_();
  // build system matrix for the FEM laplace equation
  void build_jacobian_();
  // fill out element_numbering and node_numbering
  void numberNodesEndElements_(std::vector<int> &element_types,
                        std::vector<std::vector<std::size_t> > & element_tags,
                        const std::vector<std::vector<std::size_t> > &node_tags);
  // build vectors for shape functions and fill them with initial guess values
  void initial_guess_();

 private:
  const mesh::Cell & _cell;
  std::unordered_map<size_t, size_t> _cell_numbering;
  std::unordered_map<size_t, size_t> _node_numbering;
  Eigen::SparseMatrix<double,Eigen::RowMajor> _system_matrix;
};

}  // end namespace discretization
