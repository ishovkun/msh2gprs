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
  TributaryRegion2dBase(PolyhedralElementBase & element);
  virtual ~TributaryRegion2dBase() = default;
  std::vector<std::vector<angem::Polygon<double>>> & get() {return _tributary;}
  const std::vector<std::vector<angem::Polygon<double>>> & get() const {return _tributary;}

 protected:
  PolyhedralElementBase & _element;
  // vector: face index - tributary region index - polygon
  std::vector<std::vector<angem::Polygon<double>>> _tributary;

};

}  // end namespace discretization

#endif
