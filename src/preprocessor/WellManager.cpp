#include "WellManager.hpp"
#include "angem/Collisions.hpp"
#include <numeric>  // std::accumulate
#include "VTKWriter.hpp" // debugging, provides io::VTKWriter

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
      setup_simple_well_(well);
    else
    {
      throw std::invalid_argument("Segmented wells not set up yet");
      setup_segmented_well_(well);
    }

    compute_well_index_(well);
    if (well.connected_volumes.empty())
      throw std::runtime_error("Well " + well.name + " has no conncted volumes");
    _data.wells.push_back( std::move(well) );
  }
}

void WellManager::setup_simple_well_(Well & well)
{
  const auto & grid = _data.grid;
  well.reference_depth = -well.coordinate.z();
  _well_connected_cells.emplace_back();
  std::unordered_set<size_t> processed_cells;
  // well assigned with a single coordinate
  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
    if (cell->ultimate_parent() != *cell || _edfm_method == EDFMMethod::compartmental )
    {
      const bool intersection_found = setup_simple_well_element_(well, cell->index());
      if (intersection_found)
        setup_simple_well_to_fracture_(well, cell->index());
    }
    else if (processed_cells.find(cell->ultimate_parent().index()) == processed_cells.end())
    {
      setup_simple_well_element_(well, cell->ultimate_parent().index());
    }
  if ( well.connected_volumes.empty() )
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
  well.indices.resize(well.connected_volumes.size());
  for (std::size_t i = 0; i<well.connected_volumes.size(); ++i)
  {
    const auto & cv = _data.cv_data[well.connected_volumes[i]];
    if (cv.type == discretization::ControlVolumeType::cell)  // matrix
      compute_WI_matrix_(well, i);
    else // fracture
      compute_WI_frac_(well, well.connected_volumes[i]);
  }

  const double wi_sum = std::accumulate(well.indices.begin(), well.indices.end(), 0.0);
  if (wi_sum == 0)
    throw std::runtime_error( "All well indices are zero for well: " + well.name );
}

angem::Point<3,double> WellManager::get_bounding_box_(const std::size_t icell) const
{
  const auto cell_poly = _data.grid.cell(icell).polyhedron();
  angem::Point<3,double> dir, result;
  dir = {1, 0, 0};
  result[0] = get_bounding_interval(dir, cell_poly.get());
  dir = {0, 1, 0};
  result[1] = get_bounding_interval(dir, cell_poly.get());
  dir = {0, 0, 1};
  result[2] = get_bounding_interval(dir, cell_poly.get());
  return result;
}

bool WellManager::setup_simple_well_element_(Well & well, size_t cell_index)
{
  const auto & cell = _data.grid.cell(cell_index);
  const auto p_poly_cell = cell.polyhedron();
  // define sufficiantly-large artificial segment to compute intersection with cell
  const Point cell_center = p_poly_cell->center();
  const double h = cell_center.distance(p_poly_cell->get_points()[0]);
  Point p1 = well.coordinate;
  Point p2 = well.coordinate;
  p1.z() = cell_center.z() + h * 10;
  p2.z() = cell_center.z() - h * 10;

  static const Point direction = {0, 0, -1};
  std::vector<Point> section_data;
  const double tol = 1e-6*fabs(p1.z()-p2.z());
  if (angem::collision(p1, p2, *p_poly_cell, section_data, tol))
  {
    assert( section_data.size() == 2 );

    well.connected_volumes.push_back(_dofs.cell_dof(cell_index));
    _well_connected_cells.back().push_back(cell_index);
    const double segment_length = section_data[0].distance(section_data[1]);

    well.segment_length.push_back(segment_length);
    well.directions.push_back(direction);
    // for visulatization
    _data.well_vertex_indices.emplace_back();
    _data.well_vertex_indices.back().first = _data.well_vertices.insert(section_data[0]);
    _data.well_vertex_indices.back().second = _data.well_vertices.insert(section_data[1]);
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
      if (well.force_fracture_connection())
      {
        std::cout << "fuck" << std::endl;
        abort();
      }
      else
      {
        const auto poly = face->polygon();
        std::vector<Point> section_data;
        const Point direction(0, 0, 1);
        Point p1 = well.coordinate;
        Point p2 = well.coordinate;
        p1.z() = poly.support(direction).z();
        p2.z() = poly.support(-direction).z();
        angem::collision(p1, p2, poly.plane(), section_data);
        if (!section_data.empty())
        {
          if (poly.point_inside(section_data.front()))
          {
            const size_t frac_dof = _dofs.face_dof(face->index());
            well.connected_volumes.push_back(frac_dof);
            const auto & frac_cv = _data.cv_data[frac_dof];
            std::array<double,3> bbox;
            double segment_length;
            if (section_data.size() == 1)
            {
              // find angle alpha between well direcion and its projection onto frac plane
              // segment_length = 2*r_well / sin(alpha) -- stole from Robin Hui
              const auto & plane = poly.plane();
              Point tangent1 = plane.project_vector({0, 0, 1}).normalize();
              const double cos_alpha = direction.dot(tangent1);
              segment_length = well.radius / std::sqrt( 1 - cos_alpha*cos_alpha );
              // find a bounding box for the intersection
              angem::Line<3,double> l1(section_data.front(), tangent1);
              const auto tangent2 = tangent1.cross(plane.normal()).normalize();
              angem::Line<3,double> l2(section_data.front(), tangent2);
              std::vector<Point>  sect1, sect2;
              angem::collision(l1, poly, sect1);
              angem::collision(l2, poly, sect2);
              assert( sect1.size() == 2 );
              assert( sect2.size() == 2 );
              std::array<double,3> bounding_box = {sect1.front().distance(sect1.back()),
                                                   sect2.front().distance(sect2.back()),
                                                   0};
            }
            else
            {
              segment_length = section_data.back().distance(section_data.front());
            }

            well.segment_length.push_back(segment_length);
            well.directions.push_back(direction);
            // IO::VTKWriter::write_geometry(_data.grid, cell,
            //                               "output/geometry-" + std::to_string(cell.index()) + ".vtk");
          }
        }

      }
    }
  }
 
}

void WellManager::compute_WI_matrix_(Well & well, const size_t segment)
{
  const std::size_t icell = _well_connected_cells.back()[segment];
  // const angem::Tensor2<3,double> perm = _data.get_permeability(icell);
  const auto & perm = _data.cv_data[well.connected_volumes[segment]].permeability;
  angem::Point<3,double> dx_dy_dz = get_bounding_box_(icell);
  angem::Point<3,double> productivity;
  // project on plane yz
  productivity[0] = compute_productivity(perm(1,1), perm(2,2), dx_dy_dz[1], dx_dy_dz[2],
                                         well.segment_length[segment]*fabs(well.directions[segment][0]),
                                         well.radius);
  // project on plane xz
  productivity[1] = compute_productivity(perm(0,0), perm(2,2), dx_dy_dz[0], dx_dy_dz[2],
                                         well.segment_length[segment]*
                                         fabs(well.directions[segment][1]),
                                         well.radius);
  // project on plane xy
  productivity[2] = compute_productivity(perm(0,0), perm(1,1), dx_dy_dz[0], dx_dy_dz[1],
                                         well.segment_length[segment]
                                         *fabs(well.directions[segment][2]),
                                           well.radius);
  const double wi = productivity.norm();
  if (productivity[0] < 0 or productivity[1] < 0 or productivity[2] < 0 or std::isnan(wi))
  {
    const std::string msg = "Wrong directional WI for well " + well.name + ": "
        + std::to_string( productivity[0] ) + " "
        + std::to_string(productivity[1]) + " "
        + std::to_string(productivity[2]);
    throw std::runtime_error(msg);
  }

    well.indices[segment] = wi;
}

void WellManager::compute_WI_frac_(Well & well, const size_t face_index)
{
 
}

}  // end namespace gprs_data
