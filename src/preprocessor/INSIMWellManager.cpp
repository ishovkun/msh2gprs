#include "INSIMWellManager.hpp"

namespace gprs_data {

INSIMWellManager::INSIMWellManager(std::vector<WellConfig> const & wells, mesh::Mesh const & grid,
                                   GridIntersectionSearcher & searcher)
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
    std::cout << "\nwell.perforated.size() = " << well.perforated.size() << std::endl;
    for (size_t j = 0; j < well.perforated.size(); ++j)
    {
      // find the center of the segment
      angem::Point<3,double> point;
      if ( well.coordinates.size() == 1 ) point = well.coordinates[0];
      else point = 0.5 * (well.coordinates[j] + well.coordinates[j+1]);

      // find cell that contains the point
      size_t const cell_idx = _searcher.find_cell(point);
      std::cout << "found cell_idx " << cell_idx << std::endl;

      // find closest vertex within the cell
      double mindist = std::numeric_limits<double>::max();
      size_t closest = std::numeric_limits<size_t>::max();
      for (size_t const v : _grid.cell(cell_idx).vertices()) {
        double const d = _grid.vertex(v).distance(point);
        std::cout << "testing vertex " << v  << " (" << _grid.vertex(v) << ")"
                  << " (dist = " << d
                  << std::endl;
        if ( d < mindist ) {
          mindist = d;
          closest = v;
        }
      }
      std::cout << "found " << " (" << closest << ")"
                << " closest to " << point
                << std::endl;

      _well_vertex_indices[iwell].push_back( closest );
    }
    assert( !_well_vertex_indices[iwell].empty() && "Could not find well vertices" );
  }

  for (size_t w = 0; w < _well_vertex_indices.size(); ++w)
  {
    std::cout << "well " << w << ": ";
    for (auto v : _well_vertex_indices[w])
      std::cout << v << " ";
    std::cout << std::endl;
  }


  exit(0);
}

}  // end namespace gprs_data
