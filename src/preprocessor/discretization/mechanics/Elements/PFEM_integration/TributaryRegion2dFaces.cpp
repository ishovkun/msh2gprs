#ifdef WITH_EIGEN
#include "TributaryRegion2dFaces.hpp"

namespace discretization {

TributaryRegion2dFaces::TributaryRegion2dFaces(PolyhedralElementBase & element, const size_t parent_face)
    : TributaryRegion2dBase(element, parent_face)
{
  const auto face_polygons = _element.host_cell().polyhedron()->get_face_polygons();
  build_tributary_shapes_face_(parent_face, face_polygons[parent_face]);
  mark_faces_();
}

void TributaryRegion2dFaces::build_tributary_shapes_face_(const size_t iface, const angem::Polygon<double> & face_poly)
{
  const angem::Point<3,double> center = face_poly.center();
  const std::vector<angem::Point<3,double>> & vertices = face_poly.get_points();
  std::vector<angem::Point<3,double>> result;
  for (size_t i = 0; i < vertices.size(); ++i)
  {
    size_t const v1 = i;
    size_t const v2 = (i+1) % vertices.size();
    std::vector<angem::Point<3,double>> triangle_vertices = {vertices[v1], vertices[v2], center};
    angem::Polygon<double> triangle(triangle_vertices);
    _tributary.push_back(triangle);
  }
}

void TributaryRegion2dFaces::mark_faces_()
{
  const std::vector<size_t> & face_indices = _element.get_face_domains()[_parent_face];
  const auto & grid = _element.get_grid();
  _faces.resize(_tributary.size());

  // optimization: in case of zero refinement level we have a full rule
  if (_element.get_config().subdivision_method == PolyhedralFEMSubdivision::refinement &&
      _element.get_config().order == 0) {
    for (size_t f = 0; f < face_indices.size(); ++f)
      _faces[f].push_back(face_indices[f]);
  }
  else {  // generic case with geometric check
    for (const size_t iface : face_indices)
    {
      const mesh::Face & face = grid.face(iface);
      const auto c = face.center();
      for (size_t region=0; region < _tributary.size(); ++region)  // tributary regions
        if (_tributary[region].point_inside(c))
          _faces[region].push_back(iface);
    }
  }

  _faces_center.resize( face_indices.size() );
  std::copy( face_indices.begin(), face_indices.end(), _faces_center.begin() );
}

}  // end namespace discretization

#endif
