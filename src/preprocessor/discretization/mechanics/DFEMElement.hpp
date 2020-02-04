#pragma once

#include "mesh/Cell.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include <Eigen/Sparse>  // provides SparseMatrix
#include <Eigen/Dense>  // provides MatrixXd, VectorXd

namespace discretization {

class DFEMElement
{
 public:
  DFEMElement(const mesh::Cell & cell);

 protected:
  void build_();
  // get geometric data from gmsh
  void build_triangulation_();
  // build system matrix for the FEM laplace equation
  void build_jacobian_();
  // fill out element_numbering and node_numbering
  void numberNodesEndElements_(std::vector<int> &element_types,
                        std::vector<std::vector<std::size_t> > & element_tags,
                        const std::vector<std::vector<std::size_t> > &node_tags);
  // build vectors for shape functions and fill them with initial guess values
  void initial_guess_();
  // just for debugging
  void debug_save_shape_functions_(const std::string fname = "shape_functions.vtk");
  // run iterative msrsb process to compute shape functions
  void run_msrsb_();
  // jacobi iteration over a single fine vertex
  double jacobi_iteration_(const size_t fine_vertex);
  void enforce_partition_of_unity_(const size_t fine_vertex,
                                   std::vector<Eigen::VectorXd> & solutions);
  void enforce_zero_on_boundary_(const size_t fine_vertex,
                                 std::vector<Eigen::VectorXd> & solutions);
  void build_support_boundaries_();
  bool in_support_boundary_(const size_t fine_vertex, const size_t parent_face) const;
  bool in_global_support_boundary_(const size_t fine_vertex) const;

 private:
  const mesh::Cell & _cell;
  std::vector<int> _element_types;
  std::vector<std::vector<std::size_t> > _element_tags;
  std::vector<std::vector<std::size_t> > _element_nodes;
  std::vector<std::size_t> _node_tags;
  std::vector<angem::Point<3,double>> _node_coord;
  std::unordered_map<size_t, size_t> _cell_numbering;
  std::unordered_map<size_t, size_t> _node_numbering;
  Eigen::SparseMatrix<double,Eigen::RowMajor> _system_matrix;
  // msrsb basis function
  std::vector<Eigen::VectorXd> _basis_functions;
  // mrsrb support boundaries
  std::vector<std::unordered_set<size_t>> _support_boundaries;
};

}  // end namespace discretization
