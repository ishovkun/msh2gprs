#pragma once

#ifdef WITH_EIGEN
#include "mesh/Mesh.hpp"  // provides mesh::Mesh
#include "../PolyhedralElementBase.hpp"
#include <vector>

namespace discretization {

/**
 * Thid is a base class for PFEM tributary regions for faces.
 * Per se it id dummy.
 */
class TributaryRegion2dBase {
 public:
  TributaryRegion2dBase(PolyhedralElementBase & element, const size_t parent_face)
    : _element(element), _parent_face(parent_face)
  {}
  virtual ~TributaryRegion2dBase() = default;
  std::vector<angem::Polygon<double>> & get() {return _tributary;}
  const std::vector<angem::Polygon<double>> & get() const {return _tributary;}

  // new api
  virtual size_t size() const noexcept {return _tributary.size();}
  inline const std::vector<std::size_t> & get(const size_t region_index) const noexcept
  {
    assert( region_index < size() && "Vertex index must be less the the number of vertices");
    return _faces[region_index];
  }

  inline const std::vector<std::size_t> & get_center() const noexcept
  {
    return _faces_center;
  }

  virtual double area(const size_t region_index) const {return _tributary[region_index].area();}

  virtual double area_total() const
  {
    return std::accumulate( _tributary.begin(), _tributary.end(), 0.0,
                            [](const double curr, const angem::Polygon<double>& poly)
                            { return curr + poly.area(); });
  }

 protected:
  PolyhedralElementBase & _element;
  const size_t _parent_face;
  std::vector<angem::Polygon<double>> _tributary;
  std::vector<std::vector<size_t>> _faces;
  std::vector<size_t> _faces_center;
};

}  // end namespace discretization

#endif
