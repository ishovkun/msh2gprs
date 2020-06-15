#ifdef WITH_EIGEN
#include "TributaryRegion2dBase.hpp"

namespace discretization {

TributaryRegion2dBase::TributaryRegion2dBase(PolyhedralElementBase & element)
    : _element(element)
{}

}  // end namespace discretization


#endif
