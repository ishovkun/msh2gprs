#pragma once

#ifdef WITH_EIGEN
#include "../PolyhedralElementBase.hpp"
#include "IntegrationRule2dBase.hpp"
#include "angem/Basis.hpp"

namespace discretization {

class IntegrationRule2dFull : public IntegrationRule2dBase {
 public:
  IntegrationRule2dFull(PolyhedralElementBase & element, const size_t parent_face,
                        const angem::Basis<3, double> & basis);
  virtual ~IntegrationRule2dFull() = default;
  FiniteElementData get() const;

 protected:
  void setup_storage_(FiniteElementData & data) const;
};


}  // end namespace discretization

#endif
