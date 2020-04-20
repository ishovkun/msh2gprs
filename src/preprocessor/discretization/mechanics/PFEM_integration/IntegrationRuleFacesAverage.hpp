#pragma once

#include "mesh/Mesh.hpp"  // provides mesh::Mesh
#include "../PolyhedralElementBase.hpp"
#include <vector>

namespace discretization {

/**
 * This class implements a # faces = # q-points integration rule for Polyhedral Finite Elements.
 * The values of shape functions are arerages within the tributary regions
 * There is currently no derivative corrections
 * Bishop, A displacement-based finite element formulation for
 * general polyhedra using harmonic shape functions (2014).
 */
class IntegrationRuleFacesAverage
{
 public:
  /**
   * Constructor.
   * Computes the locations of integration points at a target as well as
   * the values of shape functions, gradients, and weights.
   */
  IntegrationRuleFacesAverage(PolyhedralElementBase & element);

 protected:
  void build_tributary_shapes_cells_();
  angem::Polyhedron<double> create_pyramid_(const std::vector<size_t> & face,
                                            const std::vector<angem::Point<3,double>> & vertices) const;
  // compute shape function values, gradients, and weights in the
  // integration points in cells
  void compute_cell_fe_quantities_();
  // compute shape function values, gradients, and weights in the
  // integration points in a given face face
  void compute_face_fe_quantities_(const size_t parent_face);
  void setup_storage_();
  void build_tributary_shapes_face_(const size_t iface, const angem::Polygon<double> & face_poly);

  PolyhedralElementBase & _element;
  std::vector<angem::Polyhedron<double>> _pyramids;                         // tributary regions
  std::vector<std::vector<angem::Polygon<double>>> _face_triangles;
};

}  // end namespace discretization
