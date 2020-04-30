#pragma once

#include "mesh/Mesh.hpp"
#include "../PolyhedralElementBase.hpp"

namespace discretization {

/* This class implements cell integration rules but the integration point locations
 * are located at the faces
 */
class IntegrationRuleFacesFractures
{
 public:
  IntegrationRuleFacesFractures(PolyhedralElementBase & element);

  protected:
  void build_tributary_shapes_face_(const size_t iface, const angem::Polygon<double> & face_poly);
  // compute shape function values, gradients, and weights in the
  // integration points in a given face but for cells
  void compute_face_fe_quantities_(const size_t parent_face);
  // do proper resizing of storage vectors
  void setup_storage_();

  PolyhedralElementBase & _element;
  std::vector<std::vector<angem::Polygon<double>>> _face_triangles;  // face tributary regions
};

}  // end namespace discretization
