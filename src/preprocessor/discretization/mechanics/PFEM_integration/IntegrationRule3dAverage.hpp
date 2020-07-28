#pragma once

#ifdef WITH_EIGEN
#include "../PolyhedralElementBase.hpp"
#include "TributaryRegion3dBase.hpp"

namespace discretization {

/**
 * This class computes FEM quantities of a Polyhedral element based
 * on averaging the corresponding values wwithin tributary regions.
 */
class IntegrationRule3dAverage {
 public:
  IntegrationRule3dAverage(PolyhedralElementBase & element, const TributaryRegion3dBase  & tributary);
  virtual ~IntegrationRule3dAverage() = default;
 protected:
  void setup_storage_(PolyhedralElementBase & element, const TributaryRegion3dBase  & tributary);
  void compute_fe_values_(const std::vector<size_t> &cells, FEPointData & dst);

  PolyhedralElementBase & _element;
};

}  // end namespace discretization

#endif
