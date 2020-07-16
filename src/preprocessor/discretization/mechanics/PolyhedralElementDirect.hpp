#pragma once

#include "PolyhedralElementBase.hpp" // provides PolyhedralElementBase
#include "../flow/DoFNumbering.hpp"         // provides DoFNumbering
#include <Eigen/Sparse>                     // provides SparseMatrix

namespace discretization {

/**
 * This class implements polyhedral element method for a single 
 * cell. It uses Gmsh (linked during compilation) to triangulate an 
 * polyhdedral element and a custom FeValues class to build a system matrix 
 * for Poisson equation to compute FE shape functions.
 */
class PolyhedralElementDirect : public PolyhedralElementBase
{
 public:
  /**
   * Constructor.
   * Build FEM discretization of the cell.
   * This function has four major control parameters.
   *
   * The parameter PolyhedralFEMSubdivision in FiniteElementConfig structure
   * controls whether to use GMsh or a custom refinement.
   *
   * The parameter "order" in FiniteElemenConfig controls the
   * level of refinement for the polyhedral element.
   * If gmsh is used, then the order is the number of verties in each edge in the underlying grid.
   * If custom refinement is used, the order is the number of times we perform the
   * refinement.
   *
   * The parameter update_face_values should be used if the user needs data for
   * face integration, e.g. Neumann boundary conditions.
   *
   * The parameter udpate_fracture_values dictates whether the user needs
   * to additionally compute cell shape function data at face integration points.
   * These data is used in ADGPRS for DFM fractures, in order to solve the contact problem.
   *
   * Input:
   * \param[in] cell : grid cell to be discretized
   * \param[in] parent_grid : grid the discretized cell belongs to
   * \param[in] config : constains information about the type and order of refinement
   * \param[in] update_face_values : update face shape functions within face quadrature points
   * \param[in] update_fracture_values : update cell shape function in face quadrature points
   */
  PolyhedralElementDirect(const mesh::Cell & cell,
                          const mesh::Mesh & parent_grid,
                          const FiniteElementConfig & config);
  // get vector of cell integration points (where cell_data is defined)
  const std::vector<angem::Point<3,double>> & get_cell_integration_points() const {return _cell_gauss_points;}
  //  purely debugging purposes
  void debug_save_boundary_face_solution(const std::string fname) const;

 protected:
  // solve problems on faces
  void build_face_boundary_conditions_();
  // build system matrix for the face poisson problem
  void build_face_system_matrix_(const size_t parent_face,
                                 Eigen::SparseMatrix<double,Eigen::RowMajor> & face_system_matrix,
                                 const std::vector<size_t> & face_indices,
                                 const DoFNumbering & vertex_dofs);
  // impose boundary conditions on a poisson system for faces (to get bc's)
  void impose_bc_on_face_system_(const size_t parent_vertex,
                                 const DoFNumbering & vertex_dofs,
                                 Eigen::SparseMatrix<double,Eigen::RowMajor> & face_system_matrix,
                                 Eigen::VectorXd & rhs,
                                 const bool impose_on_matrix = true);
  void build_edge_boundary_conditions_(const std::vector<size_t> & parent_face_vertices,
                                       const std::vector<size_t> & faces);
  // find out which vertices reside on parent edges and compute the dirichlet values for them
  void build_edge_boundary_conditions_();
  // map element_grid vertices to parent cell faces
  std::vector<std::vector<size_t>> map_vertices_to_parent_faces_();
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
  // impose BC's on the cell laplace system (only rhs)
  void impose_boundary_conditions_(Eigen::VectorXd & rhs, const size_t ipv);

  void save_face_domains_(std::string fname);

 protected:
  std::vector<std::vector<size_t>> _support_edge_vertices;     // edge vertices for each parent vertex
  std::vector<std::vector<double>> _support_edge_values;       // edge dirichlet values for each parent vertex
  std::vector<std::vector<size_t>> _support_boundary_vertices; // face vertices for each parent vertex
  std::vector<std::vector<double>> _support_boundary_values;   // face dirichlet values for each parent vertex
  Eigen::SparseMatrix<double,Eigen::RowMajor> _system_matrix;  // 3d cell system matrix with no BC's
};

}  // end namespace discretization
