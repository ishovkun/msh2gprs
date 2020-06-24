#pragma once

#ifdef WITH_EIGEN
#include "../PolyhedralElementBase.hpp"

namespace discretization {

class IntegrationRule3dMS {
 public:
  /* This integration rule assumes that cell subdivisioins are the tributary regions.
   * It takes pointwise center values. Since this almost exactly what
   * multiscale medhos do, I call it an MS rule. */
  IntegrationRule3dMS(PolyhedralElementBase & element);
  virtual ~IntegrationRule3dMS() = default;

 protected:
  // do proper resizing of storage vectors
  void setup_storage_();

  PolyhedralElementBase & _element;
};

}  // end namespace discretization

#endif
