#pragma once

#ifdef WITH_EIGEN
#include "TributaryRegion3dBase.hpp"

namespace discretization {

class TributaryRegion3dVertices : public TributaryRegion3dBase
{
 public:
  TributaryRegion3dVertices(PolyhedralElementBase & element);
  virtual ~TributaryRegion3dVertices() = default;
  // get the volume of the tributary region index
  virtual double volume(const size_t region_index) const override {return _vol_tot / _n_parent_vertices ;}
  // returns the number of tributary regions
  virtual size_t size() const noexcept override {return _cells.size();}
  virtual double volume_center() const override {return _vol_tot;}
 protected:
  void mark_cells_();

  double _vol_tot;
  size_t _n_parent_vertices;
};

}  // end namespace discretization

#endif
