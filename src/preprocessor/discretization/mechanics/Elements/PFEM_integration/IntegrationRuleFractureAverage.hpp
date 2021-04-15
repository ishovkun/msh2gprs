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
                                 const TributaryRegion2dBase & tributary_2d,
                                 const size_t parent_face,
                                 const angem::Basis<3, double> & basis);
  virtual ~IntegrationRuleFractureAverage() = default;
  FiniteElementData get() const;

  protected:
  // do proper resizing of storage vectors
  void setup_storage_(FiniteElementData & data) const;
  void compute_fe_values_(const std::vector<size_t> &faces, FEPointData &dst) const;

  PolyhedralElementBase & _element;
  const TributaryRegion2dBase & _tributary;
  const size_t _parent_face;
  const angem::Basis<3,double> _basis;
};

}  // end namespace discretization

#endif
