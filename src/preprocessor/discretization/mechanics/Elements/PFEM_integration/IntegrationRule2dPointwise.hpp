#pragma once

#ifdef WITH_EIGEN
#include "TributaryRegion2dBase.hpp"

namespace discretization {

class IntegrationRule2dPointwise
{
 public:
  IntegrationRule2dPointwise(PolyhedralElementBase & element,
                           const std::vector<std::vector<angem::Polygon<double>>> & tributary_2d,
                           const size_t parent_face);
  virtual ~IntegrationRule2dPointwise() = default;
  FiniteElementData get() const;

 protected:
  void setup_storage_(FiniteElementData & data) const;
  // compute shape function values, gradients, and weights in the
  // integration points in a given face
  PolyhedralElementBase & _element;
  const std::vector<std::vector<angem::Polygon<double>>> & _tributary_2d;
  const size_t _parent_face;
};

}  // end namespace discretization


#endif
