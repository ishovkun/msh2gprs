#ifdef WITH_EIGEN
#include "TributaryRegion2dVertices.hpp"

namespace discretization {

using Point = angem::Point<3,double>;

TributaryRegion2dVertices::TributaryRegion2dVertices(PolyhedralElementBase & element,
                                                     const size_t parent_face)
    : TributaryRegion2dBase(element, parent_face)
{
  const auto polyhedron = _element._parent_cell.polyhedron();
  const auto face_polygons = polyhedron->get_face_polygons();
  build_tributary_shapes_face_(parent_face, face_polygons[parent_face]);
  mark_faces_(polyhedron->get_faces()[parent_face]);
}

void TributaryRegion2dVertices::
build_tributary_shapes_face_(const size_t iface, const angem::Polygon<double> & face_poly)
{
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
    _tributary.emplace_back(tributary_vertices);
  }
}

void TributaryRegion2dVertices::mark_faces_(const std::vector<std::size_t> & parent_verts)
{
  const auto & grid = _element._subgrid;
  _faces.resize( _tributary.size() );
  for (const size_t iface : _element._face_domains[_parent_face])
  {
    const auto & face = grid.face(iface);
    bool found = false;
    for (size_t ipv = 0; ipv < parent_verts.size(); ++ipv)
    {
      const size_t pv = parent_verts[ipv];
      if (face.has_vertex(pv))
        _faces[ipv].push_back(iface);
    }

    if (!found)
    {
      const auto c = face.center();
      auto it_min = std::min_element(parent_verts.begin(), parent_verts.end(),
                                     [grid, c](size_t i, size_t j) {
                                       return c.distance(grid.vertex(i)) < c.distance(grid.vertex(j));
                                     });
      const size_t ipv = std::distance(parent_verts.begin(), it_min);
      _faces[ipv].push_back(iface);
    }
  }

  const auto & faces = _element._face_domains[_parent_face];
  _faces_center.reserve( faces.size() );
  std::copy( faces.begin(), faces.end(),  std::back_inserter(_faces_center) );
}

}  // end namespace discretization


#endif
