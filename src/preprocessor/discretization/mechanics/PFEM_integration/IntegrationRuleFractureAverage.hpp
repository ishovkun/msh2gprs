#pragma once

#ifdef WITH_EIGEN
#include "../PolyhedralElementBase.hpp"

namespace discretization {

/* This class implements cell integration rules but the integration point locations
 * are located at the faces
 */
class IntegrationRuleFractureAverage
{
 public:
  IntegrationRuleFractureAverage(PolyhedralElementBase & element,
                               const std::vector<std::vector<angem::Polygon<double>>> & tributary_2d,
                               const size_t parent_face);
  virtual ~IntegrationRuleFractureAverage() = default;
  FiniteElementData get() const;

  protected:
  // do proper resizing of storage vectors
  void setup_storage_(FiniteElementData & data) const;

  PolyhedralElementBase & _element;
  const std::vector<std::vector<angem::Polygon<double>>> & _tributary_2d;
  const size_t _parent_face;
};

}  // end namespace discretization

#endif
