#include "INSIMWellManager.hpp"
#include "WellManager.hpp"  // WellManager::compute_WI_matrix
#include "logger/Logger.hpp"
#include "mesh/io/VTKWriter.hpp"

namespace gprs_data {

INSIMWellManager::INSIMWellManager(std::vector<WellConfig> const    & wells,
                                   mesh::Mesh const                 & grid,
                                   std::vector<std::vector<size_t>> & well_vertices,
                                   GridIntersectionSearcher         & searcher)
    : _config(wells)
    , _grid(grid)
    , _well_vertex_indices(well_vertices)
    , _searcher(searcher)
{
  setup_wells_();
}

void INSIMWellManager::compute_well_indices(std::vector<discretization::ControlVolumeData> const & cvs)
{
  for (auto & well : _wells)
    for (auto  & s : well.perforations)
      WellManager::compute_WI_matrix(well, s, cvs[s.dof].permeability);
}

void INSIMWellManager::setup_wells_()
{
  bool const vertex_indices_specified = !_well_vertex_indices.empty();
  if ( !vertex_indices_specified )
    _well_vertex_indices.resize( _config.size() );

  for (size_t iwell = 0; iwell < _config.size(); ++iwell)
  {
    auto const & conf = _config[iwell];
    Well well(conf);
    logging::log() << "Setting up well: " << well.name << std::endl;

    // std::vector<size_t> vertices;
    if ( !vertex_indices_specified )
      _well_vertex_indices[iwell] = find_well_vertices_(well);

    if ( _well_vertex_indices[iwell].empty() )
      throw std::runtime_error("Well " + well.name + " has no conncted volumes");

    well.reference_depth = compute_reference_depth_(well);
    create_well_perforations_( well, _well_vertex_indices[iwell] );

    _wells.push_back( std::move(well) );
  }
}

std::vector<size_t> INSIMWellManager::find_well_vertices_(Well const & well)
{
  std::vector<size_t> ans;
  for (size_t j = 0; j < well.segment_perforated.size(); ++j)
    if ( well.segment_perforated[j] )
    {
      // find the center of the segment
      angem::Point<3, double> point;
      if (well.simple() == 1)
        point = well.coordinate;
      else
        point = 0.5 * (well.segments[j].first + well.segments[j].second);

      // find cell that contains the point
      size_t const cell_idx = _searcher.find_cell(point);
      if ( cell_idx == _grid.n_cells_total() ) {
        mesh::IO::VTKWriter::write_geometry(_grid, "bad_grid.vtk");
        throw std::invalid_argument("could not find matching grid cell");
      }

      // find closest vertex within the cell
      double mindist = std::numeric_limits<double>::max();
      size_t closest = std::numeric_limits<size_t>::max();
      for (size_t const v : _grid.cell(cell_idx).vertices()) {
        double const d = _grid.vertex(v).distance(point);
        if (d < mindist) {
          mindist = d;
          closest = v;
        }
      }
      ans.push_back(closest);
    }

  return ans;
}

void INSIMWellManager::assign_dofs(discretization::DoFNumbering const & dofs)
{
  for (auto & well : _wells)
    for (auto & segment : well.perforations)
      segment.dof = dofs.vertex_dof( segment.element_id );
}

void INSIMWellManager::create_well_perforations_(Well & well, std::vector<size_t> const & well_vertices)
{
  well.perforations.resize( well_vertices.size() );

  if ( well.simple() ) {
    auto & s = well.perforations.back();
    s.element_id = well_vertices[0];
    compute_bounding_box_(s);
    // well.segments.push_back({ p1, p2 });
    s.length = s.bounding_box[2];
    s.direction = {0,0,1};
    s.perforated = true;
  }
  else {
    for (size_t i = 0, vertex = 0; i < well.segment_perforated.size(); ++i)
      if (well.segment_perforated[i])
      {
        auto & s = well.perforations[vertex];
        s.element_id = well_vertices[vertex];
        s.length = well.segments[i].first.distance( well.segments[i].second );
        compute_bounding_box_(s);
        s.direction = well.segments[i].first - well.segments[i].second;
        s.perforated = well.segment_perforated[i];
        vertex++;
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
    for (auto const & s : _wells[i].perforations)
      ans[i].push_back( s.element_id );
  return ans;
}

WellVTKGrid INSIMWellManager::get_well_vtk_data() const
{
  WellVTKGrid viz;

  for (auto const & well : _wells) {
    if ( well.simple() ) {
      /*
      ** The idea is to simply take two vertical-ish edges attached to the well vertex and
      ** form a line segment based on their coordinates.
      ** This is stupid but should suffice for the visualization purposes.
       */
      angem::Point<3,double> well_vertex_1 = well.coordinate;
      angem::Point<3,double> well_vertex_2 = well.coordinate;
      auto const attached_cells = _grid.vertex_cells( well.perforations[0].element_id );
      for (auto const * const cell : attached_cells)
        for (size_t const v : cell->vertices()) {
          auto const &coord = _grid.vertex(v);
          well_vertex_1[2] = std::max( well_vertex_1.z(), coord.z() );
          well_vertex_2[2] = std::min( well_vertex_2.z(), coord.z() );
        }
      viz.indices.emplace_back();
      auto & segment = viz.indices.back();
      segment.first = viz.vertices.insert(well_vertex_1);
      segment.second = viz.vertices.insert(well_vertex_2);
    }
    else {
      // simply take line segments from user-provided data
      for (auto const & s : well.segments) {
        viz.indices.emplace_back();
        auto & segment = viz.indices.back();
        segment.first = viz.vertices.insert(s.first);
        segment.second = viz.vertices.insert(s.second);
      }
    }
  }

  return viz;
}

double INSIMWellManager::compute_reference_depth_(Well const & well) const
{
  if ( std::fabs(well.reference_depth) > 1e-10 )
    return well.reference_depth;

  if ( well.simple() )
    return - well.coordinate[2];
  else {
    double deepest = std::numeric_limits<double>::max();
    for (auto const & segment : well.segments) {
      deepest = std::min( segment.first[2], deepest );
      deepest = std::min( segment.second[2], deepest );
    }
    return -deepest;
  }
}

}  // end namespace gprs_data
