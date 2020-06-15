#ifdef WITH_EIGEN
#include "TributaryRegion2dVertices.hpp"

namespace discretization {

using Point = angem::Point<3,double>;

TributaryRegion2dVertices::TributaryRegion2dVertices(PolyhedralElementBase & element)
    : TributaryRegion2dBase(element)
{
  const auto face_polygons = _element._parent_cell.polyhedron()->get_face_polygons();
  const size_t n_faces = face_polygons.size();
  for (size_t iface=0; iface < n_faces; ++iface)
    build_tributary_shapes_face_(iface, face_polygons[iface]);
}

void TributaryRegion2dVertices::
build_tributary_shapes_face_(const size_t iface, const angem::Polygon<double> & face_poly)
{
  _tributary.emplace_back();
  auto & face_tributary_polygons = _tributary.back();
  const auto & points = face_poly.get_points();
  const size_t nv = points.size();
  const Point c = face_poly.center();
  for (size_t ivertex = 0; ivertex < nv; ++ivertex)
  {
    std::vector<Point> tributary_vertices(4);
    // first vertex
    tributary_vertices[0] = face_poly.get_points()[ivertex];
    // second vertex
    const size_t vertex1 = (ivertex < nv-1) ? ivertex + 1 : 0;
    tributary_vertices[1] = 0.5 * (tributary_vertices[0] + points[vertex1]);
    // third vertex
    const size_t vertex2 = (ivertex > 0) ? ivertex - 1 : nv - 1;
    tributary_vertices[2] = 0.5 * (tributary_vertices[0] + points[vertex2]);
    // face center
    tributary_vertices[3] = c;
    face_tributary_polygons.emplace_back(tributary_vertices);
  }
}
}  // end namespace discretization


#endif
