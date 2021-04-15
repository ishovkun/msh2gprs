#pragma once

#ifdef WITH_EIGEN
#include "TributaryRegion2dBase.hpp"

namespace  discretization {

/** This class divides the faces parent of the PFEM cell into a number of tributary regions
 * equal to the number edges in the face */
class TributaryRegion2dFaces : public TributaryRegion2dBase
{
 public:
  TributaryRegion2dFaces(PolyhedralElementBase & element, const size_t parent_face);
  virtual ~TributaryRegion2dFaces() = default;

 protected:
  void build_tributary_shapes_face_(const size_t iface, const angem::Polygon<double> & face_poly);
  void mark_faces_();
};

}  // end namespace  discretization

#endif
