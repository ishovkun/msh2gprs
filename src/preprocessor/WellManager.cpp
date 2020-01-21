#include "WellManager.hpp"
#include "angem/Collisions.hpp"

namespace gprs_data {

using Point = angem::Point<3,double>;

WellManager::WellManager(const std::vector<WellConfig> & config,
                         SimData & data,
                         const discretization::DoFNumbering & dofs)
    : m_config(config), m_data(data), m_dofs(dofs)
{}

void WellManager::setup()
{
  for (const auto & conf : m_config)
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

    std::cout << "computing WI" << std::endl;
    compute_well_index_(well);
    if (well.connected_volumes.empty())
      throw std::runtime_error("Well " + well.name + " has no conncted volumes");
    m_data.wells.push_back( std::move(well) );
  }
}

void WellManager::setup_simple_well_(Well & well)
{
  std::cout << "simple well " << well.name << std::endl;
  const Point direction = {0, 0, -1};
  const auto & grid = m_data.grid;
  well.reference_depth = -well.coordinate.z();
  m_well_connected_cells.emplace_back();
  // well assigned with a single coordinate
  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
  {
    const std::unique_ptr<angem::Polyhedron<double>> p_poly_cell = cell->polyhedron();
    // define sufficiantly-large artificial segment to compute intersection with cell
    const Point cell_center = p_poly_cell->center();
    const double h = cell_center.distance(p_poly_cell->get_points()[0]);
    Point p1 = well.coordinate;
    Point p2 = well.coordinate;
    p1.z() = cell_center.z() + h * 10;
    p2.z() = cell_center.z() - h * 10;

    std::vector<Point> section_data;
    if (angem::collision(p1, p2, *p_poly_cell, section_data, 1e-6))
    {
      if (cell->ultimate_parent() != *cell)
        throw std::runtime_error("well crosses refined cell. not implemented yet");
      assert( section_data.size() == 2 );

      well.connected_volumes.push_back(m_dofs.cell_dof(cell->index()));
      m_well_connected_cells.back().push_back(cell->index());
      well.segment_length.push_back(section_data[0].distance(section_data[1]));
      well.directions.push_back(direction);
      // for visulatization
      m_data.well_vertex_indices.emplace_back();
      m_data.well_vertex_indices.back().first = m_data.well_vertices.insert(section_data[0]);
      m_data.well_vertex_indices.back().second = m_data.well_vertices.insert(section_data[1]);
    }
  }
}

void WellManager::setup_segmented_well_(Well & well)
{
  // // setup well with segments
  // std::cout << "complex well " << well.name << std::endl;
  // const auto & grid = m_data.grid;
  // m_well_connected_cells.emplace_back();
  // for (std::size_t isegment = 0; isegment < well.segments.size(); ++isegment)
  // {
  //   auto segment = well.segments[isegment];
  //   std::vector<Point> section_data;
  //   for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
  //   {
  //     const auto p_poly_cell = cell->polyhedron();
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
  // pieceman radius
  const double r = 0.28*std::sqrt(std::sqrt(k2/k1)*dx1*dx1 +
                                  std::sqrt(k1/k2)*dx2*dx2) /
                   (std::pow(k2/k1, 0.25) + std::pow(k1/k2, 0.25));
  const double j_ind = 2*M_PI*std::sqrt(k1*k2)*length/(std::log(r/radius) + skin);
  assert(j_ind >= 0);
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
    assert ( m_well_connected_cells.back().size() > i );
    const std::size_t icell = m_well_connected_cells.back()[i];
    const angem::Tensor2<3,double> perm = m_data.get_permeability(icell);
    angem::Point<3,double> dx_dy_dz = get_dx_dy_dz_(icell);
    angem::Point<3,double> productivity;
    productivity[0] =
        compute_productivity(perm(1,1), perm(2,2), dx_dy_dz[1], dx_dy_dz[2],
                             well.segment_length[i]*fabs(well.directions[i][0]),
                             well.radius);
    productivity[1] =
        compute_productivity(perm(0,0), perm(2,2), dx_dy_dz[0], dx_dy_dz[2],
                             well.segment_length[i]*fabs(well.directions[i][1]),
                             well.radius);
    productivity[2] =
        compute_productivity(perm(0,0), perm(1,1), dx_dy_dz[0], dx_dy_dz[1],
                             well.segment_length[i]*fabs(well.directions[i][2]),
                             well.radius);
    well.indices[i] = productivity.norm();
  }
}

angem::Point<3,double> WellManager::get_dx_dy_dz_(const std::size_t icell) const
{
  const auto cell_poly = m_data.grid.cell(icell).polyhedron();
  angem::Point<3,double> dir, result;
  dir = {1, 0, 0};
  result[0] = get_bounding_interval(dir, cell_poly.get());
  dir = {0, 1, 0};
  result[1] = get_bounding_interval(dir, cell_poly.get());
  dir = {0, 0, 1};
  result[2] = get_bounding_interval(dir, cell_poly.get());
  return result;
}

}  // end namespace gprs_data
