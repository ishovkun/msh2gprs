#pragma once

#ifdef WITH_EIGEN
#include "../PolyhedralElementBase.hpp"

namespace discretization {

class IntegrationRule2dFull {
 public:
  IntegrationRule2dFull(PolyhedralElementBase & element, const size_t parent_face);
  virtual ~IntegrationRule2dFull() = default;
  FiniteElementData get() const;

 protected:
  void setup_storage_(FiniteElementData & data) const;

  PolyhedralElementBase & _element;
  const size_t _parent_face;
};


}  // end namespace discretization

#endif
