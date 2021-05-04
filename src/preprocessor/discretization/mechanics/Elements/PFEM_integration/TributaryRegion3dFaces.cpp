#ifdef WITH_EIGEN
#include "TributaryRegion3dFaces.hpp"
#include <algorithm>  // iota

namespace discretization {

using Point = angem::Point<3,double>;

TributaryRegion3dFaces::TributaryRegion3dFaces(PolyhedralElementBase & element)
    : TributaryRegion3dBase(element)
{
  if (element._config.subdivision_method == PolyhedralFEMSubdivision::refinement &&
      element._config.order == 0)
  {
    // optimization for level 0 subdivision
    auto const & grid = element.get_grid();
    size_t const n_regions = element.host_cell().faces().size();
    _cells.resize(n_regions);
    for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
    {
      auto const neibs = face->neighbors();
      if (face->marker() > 0 && face->neighbors().size() == 1){
        const size_t region = static_cast<size_t>(face->marker() - 1);
        _cells[region].push_back( face->neighbors().front()->index() );
      }
    }

    _cells_center.reserve( grid.n_active_cells() );
    std::transform(grid.begin_active_cells(), grid.end_active_cells(),
                   std::back_inserter(_cells_center),
                   [](const auto cell) {return cell.index();});
  }
  else  // generic function
  {
    /* Split a parent cell into tributary regions (pyramids) */
    const auto polyhedron = _element._parent_cell.polyhedron();
    std::vector<Point> vertices = polyhedron->get_points();
    vertices.push_back(_element._parent_cell.center());
    _element._cell_gauss_points.clear();
    for (const auto &face : polyhedron->get_faces()) {
      _tributary.push_back(create_pyramid_(face, vertices));
      _element._cell_gauss_points.push_back(_tributary.back().center());
    }

    mark_cells_();
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
    size_t const v1 = face[iv];
    size_t const v2 = face[(iv+1) % face.size()];
    pyramid_faces.push_back( {v1, v2, vertex_center} );
  }
  pyramid_faces.push_back( face );  // base

  angem::Polyhedron<double> pyramid(vertices, pyramid_faces);
  return pyramid;
}

void TributaryRegion3dFaces::mark_cells_()
{
  const auto & grid = _element._subgrid;
  _cells.resize( _tributary.size() );
  for( auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell  )
  {
    const Point c = cell->center();
    for (size_t region=0; region<_tributary.size(); ++region)  // tributary regions
      if (_tributary[region].point_inside(c))
        _cells[region].push_back(cell->index());
  }

  _cells_center.reserve( grid.n_active_cells() );
  std::transform(grid.begin_active_cells(), grid.end_active_cells(),
                 std::back_inserter(_cells_center),
                 [](const auto cell) {return cell.index();});
}

}  // end namespace discretization


#endif
