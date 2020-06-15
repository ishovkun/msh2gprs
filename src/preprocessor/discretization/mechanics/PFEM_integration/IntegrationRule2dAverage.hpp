#pragma once

#ifdef WITH_EIGEN
#include "../PolyhedralElementBase.hpp"
#include "TributaryRegion2dBase.hpp"

namespace discretization {

/**
 * This class computes FEM quantities in the faces of a Polyhedral element based
 * on averaging the corresponding values within face tributary regions.
 */
class IntegrationRule2dAverage {
 public:
  IntegrationRule2dAverage(PolyhedralElementBase & element, const TributaryRegion2dBase  & tributary);
  virtual ~IntegrationRule2dAverage() = default;
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
