#include "TributaryRegion3dFull.hpp"

namespace discretization {

TributaryRegion3dFull::TributaryRegion3dFull(PolyhedralElementBase & element)
    : TributaryRegion3dBase(element)
{
  auto & grid = element.get_grid();
  size_t const n_cell = grid.n_active_cells();
  _cells.resize( n_cell );
  for (auto [cell, region] = std::tuple(grid.begin_active_cells(), 0);
       cell != grid.end_active_cells(); ++cell, ++region)
  {
    _cells[region].push_back( cell->index() );
  }
}

}  // end namespace discretization
