#include "INSIMWellManager.hpp"
#include "WellManager.hpp"  // WellManager::compute_WI_matrix

namespace gprs_data {

INSIMWellManager::INSIMWellManager(std::vector<WellConfig> const & wells, mesh::Mesh const & grid,
                                   GridIntersectionSearcher & searcher)
    : _config(wells)
    , _grid(grid)
    , _searcher(searcher)
{
  setup_wells_();
}

void INSIMWellManager::compute_well_indices(std::vector<discretization::ControlVolumeData> const & cvs)
{
  for (auto & well : _wells)
    for (auto  & s : well.segment_data)
      WellManager::compute_WI_matrix(well, s, cvs[s.dof].permeability);
}

void INSIMWellManager::setup_wells_()
{
  for (const auto & conf : _config)
  {
    Well well(conf);
    std::cout << "setting up " << well.name << std::endl;
    auto const well_vertices = find_well_vertices_(well);
    if ( well_vertices.empty() )
      throw std::runtime_error("Well " + well.name + " has no conncted volumes");

    create_well_perforations_(well, well_vertices);
  }
}

std::vector<size_t> INSIMWellManager::find_well_vertices_(Well const & well)
{
  std::vector<size_t> ans;
  for (size_t j = 0; j < well.perforated.size(); ++j)
  {
      // find the center of the segment
      angem::Point<3,double> point;
      if ( well.simple() == 1 ) point = well.coordinate;
      else point = 0.5 * (well.segments[j].first + well.segments[j].second);

      // find cell that contains the point
      size_t const cell_idx = _searcher.find_cell(point);

      // find closest vertex within the cell
      double mindist = std::numeric_limits<double>::max();
      size_t closest = std::numeric_limits<size_t>::max();
      for (size_t const v : _grid.cell(cell_idx).vertices())
      {
        double const d = _grid.vertex(v).distance(point);
        if ( d < mindist ) {
          mindist = d;
          closest = v;
        }
      }
      ans.push_back( closest );
    }

  return ans;
}

void INSIMWellManager::find_well_vertices_()
{
  _well_vertex_indices.resize( _config.size() );

  for (size_t iwell = 0; iwell < _config.size(); ++iwell)
  {
    auto const & well = _config[iwell];
    for (size_t j = 0; j < well.perforated.size(); ++j)
    {
      // find the center of the segment
      angem::Point<3,double> point;
      if ( well.coordinates.size() == 1 ) point = well.coordinates[0];
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
    assert( !_well_vertex_indices[iwell].empty() && "Could not find well vertices" );
  }
}

void INSIMWellManager::assign_dofs(discretization::DoFNumbering const & dofs)
{
  for (auto & well : _wells)
    for (auto & segment : well.segment_data)
      segment.dof = dofs.vertex_dof( segment.element_id );
}

void INSIMWellManager::create_well_perforations_(Well & well, std::vector<size_t> const & well_vertices)
{
  well.segment_data.resize( well_vertices.size() );

  if ( well.simple() ) {
    auto & s = well.segment_data.front();
    s.element_id = well_vertices[0];
    compute_bounding_box_(s);
    s.length = s.bounding_box[2];
    s.direction = {0,0,1};
  }
  else {
    for (size_t i = 0; i < well_vertices.size(); ++i) {
      auto & s = well.segment_data[i];

      if ( !well.perforated[i] )
        throw std::runtime_error("write code to handle segments without perforations");

      s.element_id = well_vertices[i];
      s.length = well.segments[i].first.distance( well.segments[i].second );
      compute_bounding_box_(s);
      s.direction = well.segments[i].first - well.segments[i].second;
    }
  }
}

void INSIMWellManager::compute_bounding_box_(discretization::WellSegment & s)
{
  auto const attached_cells = _grid.vertex_cells( s.element_id );
  double const upper = std::numeric_limits<double>::max();
  double const lower = std::numeric_limits<double>::lowest();
  angem::Point<3,double> bbox_min = {upper, upper, upper};
  angem::Point<3,double> bbox_max = {lower, lower, lower};
  for (auto const * const cell : attached_cells)
    for (size_t const v : cell->vertices()) {
      auto const &coord = _grid.vertex(v);
      for (size_t i = 0; i < 3; ++i) {
        bbox_min[i] = std::min( bbox_min[i], coord[i] );
        bbox_max[i] = std::max( bbox_max[i], coord[i] );
      }
    }

  for (size_t i = 0; i < 3; ++i)
    s.bounding_box[i] = bbox_max[i] - bbox_min[i];
}

std::vector<std::vector<size_t>> INSIMWellManager::get_well_vertices() const
{
  std::vector<std::vector<size_t>> ans( _wells.size() );
  for (size_t i = 0; i < _wells.size(); ++i)
    for (auto const & s : _wells[i].segment_data)
      ans[i].push_back( s.element_id );
  return ans;
}

}  // end namespace gprs_data
