#pragma once

#include "mesh/Mesh.hpp"  // provides mesh::Mesh
#include "FiniteElementData.hpp"  // provides FiniteElmeentData
#include <Eigen/Dense>    // provides MatrixXd, VectorXd
#include <vector>

namespace discretization {

/**
 * This class implements a # faces = # q-points integration rule for Polyhedral Finite Elements.
 * The values of shape functions are computed exactly at integration points.
 * There is currently no derivative corrections
 * Bishop, A displacement-based finite element formulation for
 * general polyhedra using harmonic shape functions (2014).
 */
class IntegrationRuleFacesPointswise
{
 public:
  /**
   * Constructor.
   * Computes the locations of integration points at a target as well as
   * the values of shape functions, gradients, and weights.
   * Input:
   * \param[in] parent_cell : a target cell in the global grid
   * \param[in] element_grid : a triangulation of the target cell
   * \param[in] basis_functions : nodal values of the harmonic basis functions
   * \param[out] cell_data : a container the output values at target cell
   * \param[out] face_data : a container the output values at the faces of the target cell
   */
  IntegrationRuleFacesPointswise(const mesh::Cell & parent_cell,
                                 const mesh::Mesh & element_grid,
                                 const std::vector<Eigen::VectorXd> & basis_functions,
                                 FiniteElementData & cell_data,
                                 std::vector<FiniteElementData> & face_data);

 protected:
  void find_integration_points_();
  // compute shape function values, gradients, and weights in the
  // integration points in cells
  void compute_cell_fe_quantities_();
  // compute shape function values, gradients, and weights in the
  // integration points in a given face face
  void compute_face_fe_quantities_(const size_t parent_face);

  const mesh::Cell & _parent_cell;                             // reference to the discretized cell
  const mesh::Mesh _element_grid;                                    // triangulation of the discretized cell
  const std::vector<Eigen::VectorXd> _basis_functions;               // numerical shape function values
  FiniteElementData _cell_data;                                // FEM values and gradients in cell integration points
  std::vector<FiniteElementData> _face_data;                                // FEM values and gradients in face integration points
};

}  // end namespace discretization
