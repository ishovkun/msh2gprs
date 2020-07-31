#pragma once

#ifdef WITH_EIGEN
#include "../PolyhedralElementBase.hpp"
#include "IntegrationRule2dBase.hpp"

namespace discretization {

/**
 * This class computes FEM quantities in the faces of a Polyhedral element based
 * on averaging the corresponding values within face tributary regions.
 */
class IntegrationRule2dAverage : public IntegrationRule2dBase
{
 public:
  IntegrationRule2dAverage(PolyhedralElementBase & element,
                           const TributaryRegion2dBase & tributary_2d,
                           const size_t parent_face,
                           const angem::Basis<3, double> & basis);
  virtual ~IntegrationRule2dAverage() = default;
  virtual FiniteElementData get() const override;
 
 protected:
  void setup_storage_(FiniteElementData & data) const;
  void compute_fe_values_(const std::vector<size_t> &faces, FEPointData &dst) const;

  const TributaryRegion2dBase & _tributary;
};

}  // end namespace discretization

#endif
