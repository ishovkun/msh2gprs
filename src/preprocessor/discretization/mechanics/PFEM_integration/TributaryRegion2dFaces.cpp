#ifdef WITH_EIGEN
#include "TributaryRegion2dFaces.hpp"

namespace discretization {

TributaryRegion2dFaces::TributaryRegion2dFaces(PolyhedralElementBase & element)
    : TributaryRegion2dBase(element)
{
  const auto face_polygons = _element._parent_cell.polyhedron()->get_face_polygons();
  const size_t n_faces = face_polygons.size();
  _tributary.resize(n_faces);
  for (size_t iface=0; iface < n_faces; ++iface)
    build_tributary_shapes_face_(iface, face_polygons[iface]);
}

void TributaryRegion2dFaces::build_tributary_shapes_face_(const size_t iface, const angem::Polygon<double> & face_poly)
{
  const angem::Point<3,double> center = face_poly.center();
  const std::vector<angem::Point<3,double>> & vertices = face_poly.get_points();
  std::vector<angem::Point<3,double>> result;
  for (size_t i=0; i<vertices.size(); ++i)
  {
    size_t v1 = i, v2;
    if ( i <  vertices.size() - 1) v2 = i + 1;
    else                           v2 = 0;
    std::vector<angem::Point<3,double>> triangle_vertices = {vertices[v1], vertices[v2], center};
    angem::Polygon<double> triangle(triangle_vertices);
    _tributary[iface].push_back(triangle);
  }
}
}  // end namespace discretization

#endif
