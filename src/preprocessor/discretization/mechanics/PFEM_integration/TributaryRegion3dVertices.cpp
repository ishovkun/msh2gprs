#ifdef WITH_EIGEN
#include "TributaryRegion3dVertices.hpp"
#include "../EdgeComparison.hpp"  // provides Edge and EdgeComparison
#include <algorithm>              // iota

namespace discretization {

using Point = angem::Point<3,double>;

TributaryRegion3dVertices::TributaryRegion3dVertices(PolyhedralElementBase & element)
    : TributaryRegion3dBase(element)
{
  mark_cells_();
}

void TributaryRegion3dVertices::mark_cells_()
{
  const size_t npv = _element._parent_cell.n_vertices();
  const auto & grid = _element._subgrid;
  _cells.resize( npv );
  std::vector<size_t> pv(npv);
  std::iota(pv.begin(), pv.end(), 0);
  for( auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell  )
  {
    bool found = false;
    for (const size_t v : cell->vertices())
      if (v < npv)
        _cells[v].push_back(cell->index());
    if (!found)
    {
      const auto c = cell->center();
      auto it_min = std::min_element(pv.begin(), pv.end(),
                                     [grid, c](size_t i, size_t j) {
                                       return c.distance(grid.vertex(i)) < c.distance(grid.vertex(j));
                                     });
      _cells[*it_min].push_back(cell->index());
    }
  }

  _cells_center.reserve( grid.n_active_cells() );
  std::transform(grid.begin_active_cells(), grid.end_active_cells(),
                 std::back_inserter(_cells_center),
                 [](const auto cell) {return cell.index();});

  _n_parent_vertices = npv;
  _vol_tot = _element._parent_cell.volume();
}


}  // end namespace discretization

#endif
