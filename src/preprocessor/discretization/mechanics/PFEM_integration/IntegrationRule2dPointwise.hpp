#pragma once

#ifdef WITH_EIGEN
#include "TributaryRegion2dBase.hpp"

namespace discretization {

class IntegrationRule2dPointwise
{
 public:
  IntegrationRule2dPointwise(PolyhedralElementBase & element, const TributaryRegion2dBase  & tributary);
  virtual ~IntegrationRule2dPointwise() = default;

 protected:
  void setup_storage_(PolyhedralElementBase & element, const TributaryRegion2dBase  & tributary);
  // compute shape function values, gradients, and weights in the
  // integration points in a given face
  void compute_face_fe_quantities_(const size_t parent_face);

  PolyhedralElementBase & _element;
  const TributaryRegion2dBase  & _tributary;
};

}  // end namespace discretization


#endif
