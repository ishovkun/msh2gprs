#ifdef WITH_EIGEN
#include "TributaryRegion3dFaces.hpp"

namespace discretization {

using Point = angem::Point<3,double>;

TributaryRegion3dFaces::TributaryRegion3dFaces(PolyhedralElementBase & element)
    : TributaryRegion3dBase(element)
{
  /* Split a parent cell into tributary regions (pyramids) */
  const auto polyhedron = _element._parent_cell.polyhedron();
  std::vector<Point> vertices = polyhedron->get_points();
  vertices.push_back(_element._parent_cell.center());
  _element._cell_gauss_points.clear();
  for (const auto & face : polyhedron->get_faces())
  {
    _tributary.push_back(create_pyramid_(face, vertices));
    _element._cell_gauss_points.push_back( _tributary.back().center() );
  }
}

angem::Polyhedron<double>
TributaryRegion3dFaces::create_pyramid_(const std::vector<size_t> & face,
                                        const std::vector<Point> & vertices) const
{
  const size_t vertex_center = vertices.size() - 1;  // HACK: I just pushed it to this array
  const auto c = _element._parent_cell.center();
  std::vector<std::vector<size_t>> pyramid_faces;
  for (size_t iv=0; iv<face.size(); ++iv)
  {
    size_t v1, v2;
    v1 = face[iv];
    if ( iv < face.size() - 1 )
      v2 = face[iv+1];
    else
      v2 = face[0];
    pyramid_faces.push_back( {v1, v2, vertex_center} );
  }
  pyramid_faces.push_back( face );  // base

  angem::Polyhedron<double> pyramid(vertices, pyramid_faces);
  return pyramid;
}


}  // end namespace discretization


#endif
