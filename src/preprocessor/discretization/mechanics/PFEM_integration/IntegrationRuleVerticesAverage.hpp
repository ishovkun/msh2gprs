#pragma once

#ifdef WITH_EIGEN
#include "mesh/Mesh.hpp"  // provides mesh::Mesh
#include "../PolyhedralElementBase.hpp"
#include <vector>

namespace discretization {

/**
 * This class implements a # vertices = # q-points integration rule for Polyhedral Finite Elements.
 * The values of shape functions are arerages within the tributary regions
 * There is no derivative corrections as in
 * Bishop, A displacement-based finite element formulation for
 * general polyhedra using harmonic shape functions (2014).
 * Instead, we average within the tributary region.
 */
class IntegrationRuleVerticesAverage
{
 public:
  IntegrationRuleVerticesAverage(PolyhedralElementBase & element);
  virtual ~IntegrationRuleVerticesAverage() = default;

 protected:
  void build_tributary_shapes_cells_();

  void build_tributary_cell_faces_(const std::vector<mesh::Edge> & edges,
                                   std::vector<std::vector<size_t>> & triburary_faces,
                                   angem::PointSet<3,double> & triburaty_vertices);

  void build_tributary_shapes_face_(const size_t iface, const angem::Polygon<double> & face_poly);
  void setup_storage_();
  // compute shape function values, gradients, and weights in the
  // integration points in cells
  void compute_cell_fe_quantities_();
  // compute shape function values, gradients, and weights in the
  // integration points in a given face
  void compute_face_fe_quantities_(const size_t parent_face);

  PolyhedralElementBase & _element;  // single PFEM element
  std::vector<angem::Polyhedron<double>> _tributary3d;
  std::vector<std::vector<angem::Polygon<double>>> _tributary2d;
};


}  // end namespace discretization

#endif
