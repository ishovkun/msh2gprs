#pragma once

#ifdef WITH_EIGEN
#include "mesh/Mesh.hpp"  // provides mesh::Mesh
#include "../PolyhedralElementBase.hpp"
#include <vector>

namespace discretization {

/**
 * This is base class for PFEM tributary regions.
 * Per se it is useless.
 */
class TributaryRegion3dBase {
 public:
  TributaryRegion3dBase(PolyhedralElementBase & element);
  virtual ~TributaryRegion3dBase() = default;
  inline std::vector<angem::Polyhedron<double>> & get() {return _tributary;}
  inline const std::vector<angem::Polyhedron<double>> & get() const {return _tributary;}
  // get cells that constitute a tributary region for selected parent vertex
  // for some rules (e.g. pointwise) that would be a single cell per region
  inline const std::vector<std::size_t> & get(const size_t region_index) const noexcept
  {
    assert( region_index < size() && "Vertex index must be less the the number of vertices");
    return _cells[region_index];
  }

  inline const std::vector<std::size_t> & get_center() const noexcept
  {
    return _cells_center;
  }

  // returns the number of tributary regions
  virtual size_t size() const noexcept {return _tributary.size();}
  // get the volume of the tributary region index
  virtual double volume(const size_t region_index) const {return _tributary[region_index].volume();}

  virtual double volume_center() const
  {
    return std::accumulate( _tributary.begin(), _tributary.end(), 0.0,
                            [](const double curr, const angem::Polyhedron<double>& poly)
                            { return curr + poly.volume(); });
  }

 protected:
  PolyhedralElementBase & _element;
  std::vector<angem::Polyhedron<double>> _tributary;
  // vector: parent vertex -> list of cells
  std::vector<std::vector<size_t>> _cells;
  // vector of cells that need to be evaluated to get fe values in center
  std::vector<size_t> _cells_center;
};

}  // end namespace discretization

#endif
