#pragma once

#include "mesh/Cell.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include <Eigen/Sparse>  // provides SparseMatrix
#include <Eigen/Dense>  // provides MatrixXd, VectorXd


namespace discretization {

class PolyhedralElementDirect
{
 public:
  PolyhedralElementDirect(const mesh::Cell & cell);

 protected:
  // main method to compute shape functions
  void build_();
  // solve problems on faces
  void build_face_boundary_conditions_();
  // identify child faces that belong to each face parent
  std::vector<std::vector<size_t>> create_face_domains_();
  // identify vertices the will constitute the linear system and create dof mapping
  std::vector<size_t> create_face_vertex_dofs_(const std::vector<size_t> & face_indices);
  // build system matrix for the face poisson problem
  void build_face_system_matrix_(const size_t face_index, const std::vector<size_t> & vertex_dofs);

 private:
  const mesh::Cell & _parent_cell;
  mesh::Mesh _element_grid;
};


}  // end namespace discretization
