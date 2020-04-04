#pragma once

#include "mesh/Cell.hpp"                    // provices mesh::Cell
#include "gmsh_interface/GmshInterface.hpp" // provides GmshInterface
#include "../flow/DoFNumbering.hpp"         // provides DoFNumbering
#include "FiniteElementData.hpp"            // provides FiniteElementData
#include <Eigen/Sparse>                     // provides SparseMatrix
#include <Eigen/Dense>                      // provides MatrixXd, VectorXd

namespace discretization {

/** 
 * This class implements polyhedral element method for a single 
 * cell. It uses Gmsh (linked during compilation) to triangulate an 
 * polyhdedral element and a custom FeValues class to build a system matrix 
 * for Poisson equation to compute FE shape functions.
 **/
class PolyhedralElementDirect
{
 public:
  /**
   * Constructor.
   * Build FEM discretization of the cell.
   * Input:
   * @param  {mesh::Cell} cell : grid cell to be discretized
   */
  PolyhedralElementDirect(const mesh::Cell & cell);
  // get FE data for volume integration
  const FiniteElementData & get_cell_data() const {return _cell_data;}
  // get FE data for surface integration
  const std::vector<FiniteElementData> & get_face_data() const {return _face_data;}
  // purely debugging purposes
  void debug_save_boundary_face_solution(const std::string fname) const;
  // purely debugging purposes
  void debug_save_shape_functions_(const std::string fname) const;

 protected:
  // main method to compute shape functions
  void build_();
  // solve problems on faces
  void build_face_boundary_conditions_();
  // identify child faces that belong to each face parent
  std::vector<std::vector<size_t>> create_face_domains_();
  // build system matrix for the face poisson problem
  void build_face_system_matrix_(const size_t iface,
                                 Eigen::SparseMatrix<double,Eigen::RowMajor> & face_system_matrix,
                                 const DoFNumbering & vertex_dofs);
  // build system matrix for the face poisson problem
  void build_face_system_matrix_(const size_t parent_face,
                                 Eigen::SparseMatrix<double,Eigen::RowMajor> & face_system_matrix,
                                 const std::vector<size_t> & face_indices,
                                 const DoFNumbering & vertex_dofs);
  // get the relation between gmsh vertex ids and grid vertices
  void compute_vertex_mapping_();
  // impose boundary conditions on a poisson system for faces (to get bc's)
  void impose_bc_on_face_system_(const size_t parent_vertex,
                                 const DoFNumbering & vertex_dofs,
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
  // put the face boundary problem solution in to the container for the boundary conditions of
  // the 3D domain system
  void append_face_solution_(const size_t pv, const Eigen::VectorXd & solution,
                             const DoFNumbering & vertex_numbering);
  // construct the laplace system matrix for the cell volume laplace equation
  void build_cell_system_matrix_();
  // impose BC's and solve laplace system to get shape functions
  void compute_shape_functions_();
  // impose BC's on the cell laplace system
  void impose_boundary_conditions_(Eigen::SparseMatrix<double,Eigen::RowMajor> & mat,
                                   Eigen::VectorXd & rhs, const size_t ipv);
  // identify the locations of the gauss points for the polyhedral element
  void find_integration_points_();
  // compute shape function values, gradients, and weights in the
  // integration points in cells
  void compute_cell_fe_quantities_();
  // compute shape function values, gradients, and weights in the
  // integration points in a given face face
  void compute_face_fe_quantities_(const size_t parent_face);
  // create a pyramid element from a cell face and cell center
  // and return its center
  angem::Point<3,double>
  create_pyramid_and_compute_center_(const std::vector<size_t> & face,
                                     const std::vector<angem::Point<3,double>>  & vertices) const;

  /**
   * Split a polygon into triangles, and compute their centers
   * Parameters:
   * \param[in] poly : a polyhedron(represents a face)
   * Returns:
   * vector of triangle center coordinates
   */
  std::vector<angem::Point<3,double>>
  split_into_triangles_and_compute_center_(const angem::Polygon<double> & poly);

  void debug_print_fe_values() const;

 private:
  const mesh::Cell & _parent_cell;                             // reference to the discretized cell
  mesh::Mesh _element_grid;                                    // triangulation of the discretized cell
  std::vector<size_t> _vertex_mapping;                         // map gmsh vertex to grid vertex
  std::vector<std::vector<size_t>> _face_domains;              // child face indices for each parent face
  std::vector<std::vector<size_t>> _support_edge_vertices;     // edge vertices for each parent vertex
  std::vector<std::vector<double>> _support_edge_values;       // edge dirichlet values for each parent vertex
  std::vector<std::vector<size_t>> _support_boundary_vertices; // face vertices for each parent vertex
  std::vector<std::vector<double>> _support_boundary_values;   // face dirichlet values for each parent vertex
  Eigen::SparseMatrix<double,Eigen::ColMajor> _system_matrix;  // 3d cell system matrix with no BC's
  std::vector<Eigen::VectorXd> _basis_functions;               // numerical shape function values
  std::vector<angem::Point<3,double>> _cell_gauss_points;      // FEM gauss points
  std::vector<std::vector<angem::Point<3,double>>> _face_gauss_points; // FEM face gauss points
  FiniteElementData _cell_data;                                // FEM values and gradients in cell integration points
  std::vector<FiniteElementData> _face_data;                                // FEM values and gradients in face integration points
};

}  // end namespace discretization
