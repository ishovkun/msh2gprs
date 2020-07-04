#pragma once

#ifdef WITH_EIGEN
#include "../PolyhedralElementBase.hpp"
#include "TributaryRegion2dBase.hpp"

namespace discretization {

/* This integration rule assumes that face subdivisioins are the tributary regions.
 * It takes pointwise sub-face center values.  */
class IntegrationRuleFractureFull {
 public:
  IntegrationRuleFractureFull(PolyhedralElementBase & element, const size_t parent_face);
  virtual ~IntegrationRuleFractureFull() = default;
  FiniteElementData get() const;

 protected:
  // do proper resizing of storage vectors
  void setup_storage_(FiniteElementData & data) const;

  PolyhedralElementBase & _element;
  const size_t _parent_face;
};


}  // end namespace discretization

#endif
