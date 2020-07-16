#pragma once

#ifdef WITH_EIGEN
#include "../PolyhedralElementBase.hpp"

namespace discretization {

/**
 * This class computes FEM quantities in the faces of a Polyhedral element based
 * on averaging the corresponding values within face tributary regions.
 */
class IntegrationRule2dAverage {
 public:
  IntegrationRule2dAverage(PolyhedralElementBase & element,
                           const std::vector<std::vector<angem::Polygon<double>>> & tributary_2d,
                           const size_t parent_face,
                           const angem::Basis<3, double> & basis);
  virtual ~IntegrationRule2dAverage() = default;
  FiniteElementData get() const;
 
 protected:
  void setup_storage_(FiniteElementData & data) const;

  PolyhedralElementBase & _element;
  const std::vector<std::vector<angem::Polygon<double>>> & _tributary_2d;
  const size_t _parent_face;
  const angem::Basis<3,double> _basis;
};

}  // end namespace discretization

#endif
