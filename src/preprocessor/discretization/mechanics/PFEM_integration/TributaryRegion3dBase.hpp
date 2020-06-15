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

 protected:
  PolyhedralElementBase & _element;
  std::vector<angem::Polyhedron<double>> _tributary;
};

}  // end namespace discretization

#endif
