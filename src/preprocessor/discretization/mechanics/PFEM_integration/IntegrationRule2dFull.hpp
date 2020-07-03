#pragma once

#ifdef WITH_EIGEN
#include "../PolyhedralElementBase.hpp"

namespace discretization {

class IntegrationRule2dFull {
 public:
  IntegrationRule2dFull(PolyhedralElementBase & element);
  virtual ~IntegrationRule2dFull() = default;

 protected:
  void setup_storage_();
  // compute shape function values, gradients, and weights in the
  // integration points in a given face
  void compute_face_fe_quantities_(const size_t parent_face);

  PolyhedralElementBase & _element;
};


}  // end namespace discretization

#endif
