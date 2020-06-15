#pragma once

#ifdef WITH_EIGEN
#include "../PolyhedralElementBase.hpp"
#include "TributaryRegion3dBase.hpp"

namespace discretization {

/**
 * This class computes FEM quntities of a Polyhedral element exactly in
 * the centers of the tributary regions. */
class IntegrationRule3dPointwise
{
 public:
   IntegrationRule3dPointwise(PolyhedralElementBase & element, const TributaryRegion3dBase  & tributary);
   virtual ~IntegrationRule3dPointwise() = default;
 protected:
  void setup_storage_();

  PolyhedralElementBase & _element;
  const TributaryRegion3dBase  & _tributary;
 };

}  // end namespace discretization

#endif
