#pragma once

#ifdef WITH_EIGEN
#include "../PolyhedralElementBase.hpp"
#include "TributaryRegion2dBase.hpp"

namespace discretization {

/* This integration rule assumes that face subdivisioins are the tributary regions.
 * It takes pointwise sub-face center values.  */
class IntegrationRuleFractureFull {
 public:
  IntegrationRuleFractureFull(PolyhedralElementBase & element);
  virtual ~IntegrationRuleFractureFull() = default;

 protected:
  //  integration points in a given face but for cells
  void compute_face_fe_quantities_(const size_t parent_face);
  // do proper resizing of storage vectors
  void setup_storage_();

  PolyhedralElementBase & _element;
};


}  // end namespace discretization

#endif