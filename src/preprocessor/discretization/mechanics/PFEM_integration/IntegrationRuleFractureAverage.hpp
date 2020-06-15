#pragma once

#ifdef WITH_EIGEN
#include "../PolyhedralElementBase.hpp"
#include "TributaryRegion2dBase.hpp"

namespace discretization {

/* This class implements cell integration rules but the integration point locations
 * are located at the faces
 */
class IntegrationRuleFractureAverage
{
 public:
  IntegrationRuleFractureAverage(PolyhedralElementBase & element, const TributaryRegion2dBase  & tributary);
  virtual ~IntegrationRuleFractureAverage() = default;
  protected:
  // compute shape function values, gradients, and weights in the
  // integration points in a given face but for cells
  void compute_face_fe_quantities_(const size_t parent_face);
  // do proper resizing of storage vectors
  void setup_storage_();

  PolyhedralElementBase & _element;
  const TributaryRegion2dBase  & _tributary;
};

}  // end namespace discretization

#endif
