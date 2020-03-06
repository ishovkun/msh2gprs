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
  void build_face_system_matrix_(const size_t iface,
                                 Eigen::SparseMatrix<double,Eigen::RowMajor> & face_system_matrix,
                                 const std::vector<size_t> & vertex_dofs);
  // get the relation between gmsh vertex ids and grid vertices
  void compute_vertex_mapping_();
  // impose boundary conditions on a poisson system for faces (to get bc's)
  void impose_bc_on_face_system_(const size_t parent_vertex,
                                 const std::vector<size_t> & vertex_dofs,
                                 Eigen::SparseMatrix<double,Eigen::RowMajor> & face_system_matrix,
                                 Eigen::VectorXd & rhs);
  void build_edge_boundary_conditions_(const std::vector<size_t> & parent_face_vertices,
                                       const std::vector<size_t> & faces);
  // find out which vertices reside on parent edges and compute the dirichlet values for them
  void build_edge_boundary_conditions_();
  // map element_grid vertices to parent cell faces
  std::vector<std::vector<size_t>> map_vertices_to_parent_faces_();
  // map parent cell vertices to parent cell faces
  std::vector<std::list<size_t>> map_parent_vertices_to_parent_faces_();

 private:
  const mesh::Cell & _parent_cell;
  mesh::Mesh _element_grid;
  std::vector<size_t> _vertex_mapping;  // map gmsh vertex to grid vertex
  // data for constructing face boundary-value problem
  std::vector<std::vector<size_t>> _support_edge_vertices;  // edge vertices for each parent vertex
  std::vector<std::vector<double>> _support_edge_values;    // edge dirichlet values for each parent vertex

};


}  // end namespace discretization
