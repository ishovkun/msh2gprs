#include "INSIMWellManager.hpp"

namespace gprs_data {

INSIMWellManager::INSIMWellManager(std::vector<WellConfig> const & wells, mesh::Mesh const & grid,
                                   GridIntersectionSearcher const & searcher)
    : _wells(wells)
    , _grid(grid)
    , _searcher(searcher)
{
  find_well_vertices_();
}

void INSIMWellManager::find_well_vertices_()
{
  _well_vertex_indices.resize( _wells.size() );

  for (size_t iwell = 0; iwell < _wells.size(); ++iwell)
  {
    auto const & well = _wells[iwell];
    for (size_t j = 0; j < well.perforated.size(); ++j) {
      // find the center of the segment
      angem::Point<3,double> point;
      if ( well.coordinates.size() == 1 )  // simple well
        point = well.coordinates[0];
      else point = 0.5 * (well.coordinates[j] + well.coordinates[j+1]);

      // find cell that contains the point
      size_t const cell_idx = _searcher.find_cell(point);

      // find closest vertex within the cell
      double mindist = std::numeric_limits<double>::max();
      size_t closest = std::numeric_limits<size_t>::max();
      for (size_t const v : _grid.cell(cell_idx).vertices()) {
        double const d = _grid.vertex(v).distance(point);
        if ( d < mindist ) {
          mindist = d;
          closest = v;
        }
      }

      _well_vertex_indices[iwell].push_back( closest );
    }
    assert( !_well_vertex_indices.empty() && "Could not find well vertices" );
  }
}

}  // end namespace gprs_data
