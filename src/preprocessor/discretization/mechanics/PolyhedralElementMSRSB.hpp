#pragma once

#include "PolyhedralElementBase.hpp" // provides PolyhedralElementBase
#include "JacobiPreconditioner.hpp"  // provides JacobiPreconditioner
#include <Eigen/Sparse>              // provides SparseMatrix

namespace discretization {

/**
 * This class implements polyhedral element method for a single
 * cell.
 * It uses Gmsh (linked during compilation) or a custom refinement class to triangulate an
 * polyhdedral element and a custom FeValues class to build a system matrix
 * for Poisson equation to compute FE shape functions.
 * Instead of a linear solver it uses MSRSB with the Jacobi-smoothing strategy
 */
class PolyhedralElementMSRSB : public PolyhedralElementBase
{
 public:
  /* Constructor */
  PolyhedralElementMSRSB(const mesh::Cell & cell, const mesh::Mesh & grid, const FiniteElementConfig & config);
  // get FE data for volume integration
  const FiniteElementData & get_cell_data() const;
  // get FE data for surface integration
  const FiniteElementData & get_face_data() const;
  // just for debugging
  void debug_save_shape_functions_(const std::string fname = "shape_functions.vtk");
  // get vector of cell integration points (where cell_data is defined)
  const std::vector<angem::Point<3,double>> & get_cell_integration_points() const {return _cell_gauss_points;}

 protected:
  void build_();
  // build system matrix for the FEM laplace equation
  void build_jacobian_();
  void compute_shape_functions_();
  void impose_boundary_conditions_(Eigen::SparseMatrix<double,Eigen::RowMajor> & mat,
                                   Eigen::VectorXd & rhs,
                                   const size_t ipv);
  // fill out element_numbering and node_numbering
  void numberNodesEndElements_(std::vector<int> &element_types,
                        std::vector<std::vector<std::size_t> > & element_tags,
                        const std::vector<std::vector<std::size_t> > &node_tags);
  void initialize_shape_functions_();
  // build vectors for shape functions and fill them with initial guess values
  void initial_guess_();
  // different guess
  void initial_guess2_();
  // run iterative msrsb process to compute shape functions
  void run_msrsb_();
  // jacobi iteration over a single fine vertex
  double jacobi_iteration_(std::vector<Eigen::VectorXd> & solutions,
                           const JacobiPreconditioner & prec);
  void enforce_partition_of_unity_(const size_t fine_vertex,
                                   std::vector<Eigen::VectorXd> & solutions);
  void enforce_zero_on_boundary_(const size_t fine_vertex,
                                 std::vector<Eigen::VectorXd> & solutions);
  void build_support_boundaries_();
  // mark vertices and compute the value of the dirichlet BC for it
  void build_support_edges_();
  // conveniance function to quickly check if a fine vertex is on support boundary of
  // the parent vertex
  bool in_support_boundary_(const size_t fine_vertex, const size_t parent_node) const;
  bool in_global_support_boundary_(const size_t fine_vertex) const;
  // debug function to help visualize the support regions
  void save_support_boundaries_();
  // debug function to help visualize the support regions
  void save_support_edges_();
  // identify the locations of the gauss points for the dfem element
  void find_integration_points_();
  // compute shape function values, gradients, and weights in the
  // integration points
  void compute_fe_quantities_();
  // create a pyramid element from a cell face and cell center
  // and return its center
  angem::Polyhedron<double>
  create_pyramid_(const std::vector<size_t> & face,
                  const std::vector<angem::Point<3,double>>  & vertices) const;

 private:
  Eigen::SparseMatrix<double,Eigen::RowMajor> _system_matrix;
  // mrsrb support boundaries
  std::vector<std::unordered_set<size_t>> _support_boundaries;
  std::vector<std::vector<size_t>> _support_edge_vertices;
  std::vector<std::vector<double>> _support_edge_values;
};

}  // end namespace discretization
