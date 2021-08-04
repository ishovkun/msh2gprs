#pragma once
#include "IntegrationRule3d.hpp"

namespace discretization {

class IntegrationRule3dEvaluation : public IntegrationRule3d {
 public:
  IntegrationRule3dEvaluation(PolyhedralElementBase const & element, TributaryRegion3dBase const  & tributary);
  virtual ~IntegrationRule3dEvaluation() = default;
};

}  // end namespace discretization
