#ifdef WITH_EIGEN
#include "TributaryRegion3dBase.hpp"

namespace discretization {

TributaryRegion3dBase::TributaryRegion3dBase(PolyhedralElementBase & element)
    : _element(element)
{}

}  // end namespace discretization

#endif
