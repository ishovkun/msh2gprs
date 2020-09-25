#include "WellManager.hpp"
#include "angem/Collisions.hpp"
#include <numeric>  // std::accumulate
#include "logger/Logger.hpp"

namespace gprs_data {

using Point = angem::Point<3,double>;

WellManager::WellManager(const std::vector<WellConfig> & config,
                         SimData & data,
                         const discretization::DoFNumbering & dofs,
                         const EDFMMethod edfm_method)
    : _config(config), _data(data), _dofs(dofs), _edfm_method(edfm_method)
{}

void WellManager::setup()
{
  for (const auto & conf : _config)
  {
    Well well(conf);
    std::cout << "setting up " << well.name << std::endl;
    if (well.simple())
    {
      // setup_simple_well_(well);
      setup_simple_well_fast_(well);
    }
    else
    {
      throw std::invalid_argument("Segmented wells not set up yet");
      setup_segmented_well_(well);
    }

    compute_well_index_(well);
    if (well.segment_data.empty())
      throw std::runtime_error("Well " + well.name + " has no conncted volumes");
    _data.wells.push_back( std::move(well) );
  }
}

void WellManager::setup_simple_well_fast_(Well & well)
{
  if (!_data.grid_searcher)
    throw std::runtime_error("Method is not supported without a valid grid searcher");
  auto & searcher = *_data.grid_searcher;
  auto p1 = well.coordinate, p2 = well.coordinate;
  well.reference_depth = -well.coordinate.z();
  p1[2] = searcher.top();
  p2[2] = searcher.bottom();
  angem::LineSegment<double> segment(p1, p2);
  std::unordered_set<size_t> processed_cells;
  for (const size_t cell_index : searcher.collision(segment))
  {
    const auto & cell = _data.grid.cell(cell_index);
    const size_t parent_index = cell.ultimate_parent().index();
    if (parent_index == cell_index || _edfm_method == EDFMMethod::compartmental )
    {
      const bool intersection_found = setup_simple_well_matrix_(well, cell_index);
      if (intersection_found)
      {
        setup_simple_well_to_fracture_(well, cell_index);
      }
    }
    else  // (p)EDFM case
    {
      // we should process merged cells only once
      //  here we pass ultimate_parent since split cells are merged back together
      if (processed_cells.find(parent_index) == processed_cells.end())
        setup_simple_well_matrix_(well, parent_index);
      // single edfm segments are not merged, we process all of them
      setup_simple_well_to_fracture_(well, cell_index);
    }
  }
  if ( well.segment_data.empty() )
    throw std::invalid_argument("Well " + well.name + " has no connected volumes");
}

void WellManager::setup_simple_well_(Well & well)
{
  const auto & grid = _data.grid;
  well.reference_depth = -well.coordinate.z();
  // _well_connected_cells.emplace_back();
  std::unordered_set<size_t> processed_cells;
  // well assigned with a single coordinate
  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
    if (cell->ultimate_parent() == *cell || _edfm_method == EDFMMethod::compartmental )
    {
      const bool intersection_found = setup_simple_well_matrix_(well, cell->index());
      if (intersection_found)
      {
        setup_simple_well_to_fracture_(well, cell->index());
      }
    }
    else  // (p)EDFM case
    {
      // we should process merged cells only once
      //  here we pass ultimate_parent since split cells are merged back together
      if (processed_cells.find(cell->ultimate_parent().index()) == processed_cells.end())
        setup_simple_well_matrix_(well, cell->ultimate_parent().index());
      // single edfm segments are not merged, we process all of them
      setup_simple_well_to_fracture_(well, cell->index());
    }
  if ( well.segment_data.empty() )
    throw std::invalid_argument("Well " + well.name + " has no connected volumes");
}

void WellManager::setup_segmented_well_(Well & well)
{
  // setup well with segments
  // std::cout << "complex well " << well.name << std::endl;
  // const auto & grid = m_data.grid;
  // m_well_connected_cells.emplace_back();
  // for (std::size_t isegment = 0; isegment < well.segments.size(); ++isegment)
  // {
  //   auto segment = well.segments[isegment];
  //   std::vector<Point> section_data;
  //   for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
  //   {
  //     const auto p_poly_cell = cell>polyhedron();
  //     if (angem::collision(segment.first, segment.second,
  //                          *p_poly_cell, section_data, 1e-6))
  //     {
  //       if (section_data.size() != 2)
  //       {
  //         std::cout << "just touching cell " << cell->index() << std::endl;
  //         section_data.clear();
  //         continue;
  //       }
  //       std::cout << "fully occupying cell " << cell->index() << std::endl;

  //       well.connected_volumes.push_back(m_data.cell_cv_indices[cell->index()]);
  //       m_well_connected_cells.back().push_back(cell->index());
  //       well.segment_length.push_back(section_data[0].distance(section_data[1]));
  //       well.directions.push_back(segment.second - segment.first);
  //       well.directions.back().normalize();

  //       // for visulatization
  //       m_data.well_vertex_indices.emplace_back();
  //       m_data.well_vertex_indices.back().first = m_data.well_vertices.insert(section_data[0]);
  //       m_data.well_vertex_indices.back().second = m_data.well_vertices.insert(section_data[1]);

  //       // auto-detect reference depth for bhp
  //       if (!well.reference_depth_set)
  //         if(cell->center()[2] < well.reference_depth)
  //           well.reference_depth = cell->center()[2];
  //       section_data.clear();
  //     }
  //   }
  // }

  // well.reference_depth_set = true;

  // // error if no connected volumes
  // if (well.connected_volumes.empty())
  //   throw std::invalid_argument("well " + well.name + " outside of the domain. aborting");
}

double compute_productivity(const double k1, const double k2,
                            const double dx1, const double dx2,
                            const double length, const double radius,
                            const double skin = 0)
{
  if (k1 == 0 or k2 == 0)
  {
    if (k1 != k2)
      throw std::runtime_error("Zero perm only in one direction is not supported");
    return 0.0;
  }
  // pieceman radius
  const double r = 0.28*std::sqrt(std::sqrt(k2/k1)*dx1*dx1 +
                                  std::sqrt(k1/k2)*dx2*dx2) /
                   (std::pow(k2/k1, 0.25) + std::pow(k1/k2, 0.25));
  double j_ind = 2*M_PI*std::sqrt(k1*k2)*length/(std::log(r/radius) + skin);
  // treatment for almost zero jind
  if (j_ind > - k1 * 1e-8 && j_ind < 0) j_ind = 0;
  // else error is thrown in the parent method
  return j_ind;
}

inline
double get_bounding_interval(const angem::Point<3,double>     & direction,
                             const angem::Polyhedron<double> * poly)
{
  assert(fabs(direction.norm() - 1) < 1e-8);
  angem::Point<3,double> neg_direction = - direction;
  return fabs((poly->support(direction) - poly->support(neg_direction)).dot(direction));
}

void WellManager::compute_well_index_(Well &well)
{
  for (auto & segment : well.segment_data)
  {
    const auto & cv = _data.cv_data[segment.dof];
    if (cv.type == discretization::ControlVolumeType::cell)  // matrix
      compute_WI_matrix_(well, segment);
    else // fracture
      compute_WI_frac_(well, segment);
  }

  // const double wi_sum = std::accumulate(well.indices.begin(), well.indices.end(), 0.0);
  const double wi_sum = std::accumulate(well.segment_data.begin(), well.segment_data.end(),
                                        0.0, [](const double sofar, const auto & segment)
                                             {return segment.wi + sofar;});
  if (wi_sum == 0)
    throw std::runtime_error( "All well indices are zero for well: " + well.name );
}

std::array<double,3> WellManager::get_bounding_box_(const std::size_t icell) const
{
  const auto cell_poly = _data.grid.cell(icell).polyhedron();
  angem::Point<3,double> dir;
  std::array<double,3> result;
  dir = {1, 0, 0};
  result[0] = get_bounding_interval(dir, cell_poly.get());
  dir = {0, 1, 0};
  result[1] = get_bounding_interval(dir, cell_poly.get());
  dir = {0, 0, 1};
  result[2] = get_bounding_interval(dir, cell_poly.get());
  return result;
}

bool WellManager::setup_simple_well_matrix_(Well & well, size_t cell_index)
{
  const auto & cell = _data.grid.cell(cell_index);
  const auto p_poly_cell = cell.polyhedron();
  Point p1 = well.coordinate;
  Point p2 = well.coordinate;
  static const Point direction = {0, 0, -1};
  const double h = cell.vertex_coordinates()[1].distance(cell.vertex_coordinates()[0]);
  p1.z() = p_poly_cell->support(direction) .z() - h;
  p2.z() = p_poly_cell->support(-direction).z() + h;

  std::vector<Point> section_data;
  const double tol = 1e-6*fabs(p1.z()-p2.z());
  if (angem::collision(angem::LineSegment<double>(p1, p2), *p_poly_cell, section_data, tol))
  {
    assert( section_data.size() == 2 );
    const double segment_length = section_data[0].distance(section_data[1]);

    // for visulatization
    _data.well_vertex_indices.emplace_back();
    _data.well_vertex_indices.back().first = _data.well_vertices.insert(section_data[0]);
    _data.well_vertex_indices.back().second = _data.well_vertices.insert(section_data[1]);

    well.segment_data.emplace_back();
    discretization::WellSegment & segment = well.segment_data.back();
    segment.dof = _dofs.cell_dof(cell_index);
    segment.element_id = cell_index;
    segment.length = segment_length;
    segment.direction = direction;
    segment.perforated = true;
    segment.bounding_box = get_bounding_box_(cell_index);
    return true;
  }
  return false;
}

void WellManager::setup_simple_well_to_fracture_(Well & well, size_t cell_index)
{
  const auto & cell = _data.grid.cell(cell_index);
  for (const auto * face : cell.faces())
  {
    if (_dofs.is_active_face(face->index()))  // it's a fracture
    {
      const auto poly = face->polygon();
      std::vector<Point> section_data;
      const Point direction(0, 0, 1);
      Point p1 = well.coordinate;
      Point p2 = well.coordinate;
      const double h = face->vertex_coordinates()[1].distance(face->vertex_coordinates()[0]);
      p1.z() = poly.support(direction).z() + h;
      p2.z() = poly.support(-direction).z() - h;
      angem::collision(p1, p2, poly.plane(), section_data);
      if (section_data.empty())
      {
        if (well.force_fracture_connection())
        {
          logging::debug() << "forcing " << well.name << " to connect to frac face " << face->index()
                           << std::endl;
          const auto center = face->center();
          p1.x() = center.x();
          p2.x() = center.x();
          p1.y() = center.y();
          p2.y() = center.y();
          angem::collision(p1, p2, poly.plane(), section_data);
          if (section_data.empty())
            throw std::runtime_error("could not connect well " + well.name +
                                     " in cell " + std::to_string(cell_index));
        }
      }
      if (!section_data.empty())
      {
        if (poly.point_inside(section_data.front(), /*tol=*/ 0.01 * h) ||
            section_data.size() == 2)
        {
          well.segment_data.emplace_back();
          auto & segment = well.segment_data.back();
          segment.dof = _dofs.face_dof(face->index());
          const auto & plane = poly.plane();
          const Point tangent1 = plane.project_vector({0, 0, 1}).normalize();
          const Point tangent2 = tangent1.cross(plane.normal()).normalize();
          // well.connected_volumes.push_back(frac_dof);
          const auto & frac_cv = _data.cv_data[segment.dof];
          if (section_data.size() == 1)
          {
            // find angle alpha between well direcion and its projection onto frac plane
            // segment_length = 2*r_well / sin(alpha) -- stole from Robin Hui
            const double cos_alpha = direction.dot(tangent1);
            segment.length = well.radius / std::sqrt( 1 - cos_alpha*cos_alpha );
            // find a bounding box for the intersection
            angem::Line<3,double> l1(section_data.front(), tangent1);
            angem::Line<3,double> l2(section_data.front(), tangent2);
            std::vector<Point>  sect1, sect2;
            angem::collision(l1, poly, sect1);
            angem::collision(l2, poly, sect2);
            assert( sect1.size() == 2 );
            assert( sect2.size() == 2 );
            segment.bounding_box = {sect1.front().distance(sect1.back()),
                                    sect2.front().distance(sect2.back()), 0};
          }
          else  // well is colinear to the segment size() == 2
          {
            segment.length = section_data.back().distance(section_data.front());
            angem::Line<3,double> l2(0.5 * (section_data.front() + section_data.back()),
                                     tangent2);
            std::vector<Point> sect2;
            angem::collision(l2, poly, sect2);
            segment.bounding_box = {segment.length, sect2.front().distance(sect2.back()), 0};
          }

          segment.direction = direction;
        }
      }
    }
  }
 
}

void WellManager::compute_WI_matrix_(Well & well, discretization::WellSegment & segment)
{
  // const std::size_t icell = _well_connected_cells.back()[segment];
  // const angem::Tensor2<3,double> perm = _data.get_permeability(icell);
  const auto & perm = _data.cv_data[segment.dof].permeability;
  const auto & dx_dy_dz = segment.bounding_box;
  angem::Point<3,double> productivity;
  // project on plane yz
  productivity[0] = compute_productivity(perm(1,1), perm(2,2), dx_dy_dz[1], dx_dy_dz[2],
                                         segment.length*fabs(segment.direction[0]),
                                         well.radius);
  // project on plane xz
  productivity[1] = compute_productivity(perm(0,0), perm(2,2), dx_dy_dz[0], dx_dy_dz[2],
                                         segment.length* fabs(segment.direction[1]),
                                         well.radius);
  // project on plane xy
  productivity[2] = compute_productivity(perm(0,0), perm(1,1), dx_dy_dz[0], dx_dy_dz[1],
                                         segment.length*fabs(segment.direction[2]),
                                         well.radius);
  segment.wi = productivity.norm();
  if (productivity[0] < 0 or productivity[1] < 0 or productivity[2] < 0 or std::isnan(segment.wi))
  {
    const std::string msg = "Wrong directional WI for well " + well.name + ": "
        + std::to_string( productivity[0] ) + " "
        + std::to_string(productivity[1]) + " "
        + std::to_string(productivity[2]);
    throw std::runtime_error(msg);
  }
}

void WellManager::compute_WI_frac_(Well & well, discretization::WellSegment & segment)
{
  // refer to the slide: fracture wellbore index of my cedfm presentation
  // fracture-well intersection is decomposed into a radial flow part and
  // a reduced linear flow part

  const auto & frac = _data.cv_data[segment.dof];
  const double perm = frac.permeability(0,0);
  // full linear slot
  const double wi_lin_slot = 8 * frac.aperture * segment.bounding_box[0] * perm;
  // reduced linear slot
  const double wi_red_slot = 8 * frac.aperture * (segment.length - 2 * well.radius) * perm;

  const double wi_rad_slot = compute_productivity(perm, perm,
                                                  segment.bounding_box[0],
                                                  segment.bounding_box[1],
                                                  frac.aperture, well.radius,
                                                  /*skin =*/ 0);
  segment.wi = std::min(wi_red_slot + wi_rad_slot, wi_lin_slot);
}

}  // end namespace gprs_data
