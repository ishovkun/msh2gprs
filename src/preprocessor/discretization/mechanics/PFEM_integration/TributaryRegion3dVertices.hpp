#pragma once

#ifdef WITH_EIGEN
#include "TributaryRegion3dBase.hpp"

namespace discretization {

class TributaryRegion3dVertices : public TributaryRegion3dBase
{
 public:
  TributaryRegion3dVertices(PolyhedralElementBase & element);
  virtual ~TributaryRegion3dVertices() = default;
 protected:
  void build_tributary_cell_faces_(const std::vector<mesh::Edge> & edges,
                                   const mesh::Face & face,
                                   std::vector<std::vector<size_t>> & tributary_faces,
                                   angem::PointSet<3,double> & tributary_vertices);
};

}  // end namespace discretization

#endif