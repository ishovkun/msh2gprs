#pragma once

#ifdef WITH_EIGEN
#include "../PolyhedralElementBase.hpp"

namespace discretization {

class IntegrationRule3dFull {
 public:
  /* This integration rule assumes that cell subdivisioins are the tributary regions.
   * It takes pointwise center values. Since this almost exactly what
   * multiscale medhos do, I call it an MS rule. */
  IntegrationRule3dFull(PolyhedralElementBase & element);
  virtual ~IntegrationRule3dFull() = default;

 protected:
  // do proper resizing of storage vectors
  void setup_storage_();

  PolyhedralElementBase & _element;
};

}  // end namespace discretization

#endif
