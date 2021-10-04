#pragma once
#ifdef WITH_EIGEN
#include "TributaryRegion3dBase.hpp"

namespace discretization {

class TributaryRegion3dFull : public TributaryRegion3dBase {
 public:
  TributaryRegion3dFull(PolyhedralElementBase & element);
  virtual ~TributaryRegion3dFull() = default;

};


}  // end namespace discretization

#endif
