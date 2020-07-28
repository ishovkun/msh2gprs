#pragma once

#ifdef WITH_EIGEN
#include "TributaryRegion3dBase.hpp"

namespace discretization {

class TributaryRegion3dFaces : public TributaryRegion3dBase
{
 public:
  TributaryRegion3dFaces(PolyhedralElementBase & element);
  virtual ~TributaryRegion3dFaces() = default;

 protected:
  angem::Polyhedron<double> create_pyramid_(const std::vector<size_t> & face,
                                            const std::vector<angem::Point<3,double>> & vertices) const;
  void mark_cells_();
};

}  // end namespace discretization

#endif
