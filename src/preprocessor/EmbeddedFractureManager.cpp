#include "EmbeddedFractureManager.hpp"
#include "angem/CollisionGJK.hpp"  // collisionGJK
#include "angem/Collisions.hpp"    // angem::split
#include <stdexcept>
#include <utility>                 // provides std::pair
#include <fstream>                 // debug
#include "VTKWriter.hpp"


namespace gprs_data {

using std::vector;
using Point = angem::Point<3,double>;
using std::pair;

EmbeddedFractureManager::EmbeddedFractureManager(
    std::vector<EmbeddedFractureConfig> &config, const EDFMMethod edfm_method,
    const double min_dist_to_node, SimData &data)
    : config(config), m_method(edfm_method),
      _min_dist_to_node(min_dist_to_node), m_data(data), m_grid(data.grid),
      _splitter(m_grid)
{
  if (_min_dist_to_node < 1e-10)
    throw std::invalid_argument("EDFM min distance to node is too small");
}

void EmbeddedFractureManager::split_cells()
{
  int face_marker = find_maximum_face_marker_() + 1;
  for (std::size_t i=0; i<config.size(); ++i)
  {
    auto & frac = config[i];
    std::vector<size_t> cells_to_split;
    // iteratively shift fracture if it collides with any grid vertices
    size_t iter = 0;
    // while (!find_edfm_cells_(*frac.body, cells_to_split))
    while (!find_edfm_cells_fast_(*frac.body, cells_to_split))
    {
      cells_to_split.clear();
      if (++iter > 100)
        throw std::runtime_error("Cannot move fracture to avoid collision with vertices");
    }
    if (cells_to_split.empty())
      throw std::invalid_argument("embedded fracture "+std::to_string(i) + " intersects zero cells");

    split_cells_(*frac.body, cells_to_split, face_marker);
    m_marker_config.insert({ face_marker, i });
    face_marker++;
  }
  // std::cout << "active cells: ";
  // int cnt = 0;
  // for (auto cell = m_grid.begin_active_cells(); cell != m_grid.end_active_cells(); ++cell)
  // {
  //   std::cout << cell->index() << "("<<cnt<<") ";
  //   if ((++cnt) % 10 == 0)
  //     std::cout << std::endl;
  // }
}

void EmbeddedFractureManager::split_cells_(angem::Polygon<double> & fracture,
                                           std::vector<size_t> & cells,
                                           const int face_marker)
{
  const auto & plane = fracture.plane();
  for (const size_t icell : cells)
  {
    mesh::Cell & old_cell = m_grid.cell(icell);
    _splitter.split_cell(old_cell, plane, face_marker);
  }
}

bool EmbeddedFractureManager::find_edfm_cells_(angem::Polygon<double> & fracture,
                                               std::vector<size_t> & cells)
{
  angem::CollisionGJK<double> collision;

  for (auto cell = m_grid.begin_active_cells(); cell != m_grid.end_active_cells(); ++cell)
  {
    const std::unique_ptr<angem::Polyhedron<double>> polyhedron = cell->polyhedron();
    if (collision.check(fracture, *polyhedron))
    {
      cells.push_back(cell->index());

      // check if some vertices are too close to the fracture
      // and if so move a fracture a little bit
      const std::vector<Point> & vertices = polyhedron->get_points();
      for (const Point & vertex : vertices)
      {
        const double h = polyhedron->center().distance(vertex);
        const double min_dist = h * _min_dist_to_node;
        if ( std::fabs( fracture.plane().signed_distance(vertex)) < min_dist )
        {
          const Point shift = h/50 * fracture.plane().normal();
          fracture.move(shift);
          return false;
        }
      }
    }
  }

  return true;
}

bool EmbeddedFractureManager::find_edfm_cells_fast_(angem::Polygon<double> & fracture, std::vector<size_t> & cells)
{
  angem::CollisionGJK<double> collision;
  for (const size_t cell_index : m_data.grid_searcher->collision(fracture))
  {
    const auto & cell = m_grid.cell(cell_index);
    const std::unique_ptr<angem::Polyhedron<double>> polyhedron = cell.polyhedron();
    if (collision.check(fracture, *polyhedron))
    {
      cells.push_back(cell.index());

      // check if some vertices are too close to the fracture
      // and if so move a fracture a little bit
      const std::vector<Point> & vertices = polyhedron->get_points();
      for (const Point & vertex : vertices)
      {
        const double h = polyhedron->center().distance(vertex);
        const double min_dist = h * _min_dist_to_node;
        if ( std::fabs( fracture.plane().signed_distance(vertex)) < min_dist )
        {
          const Point shift = h/50 * fracture.plane().normal();
          fracture.move(shift);
          return false;
        }
      }
    }
  }
  return true;
}


std::vector<DiscreteFractureConfig> EmbeddedFractureManager::generate_dfm_config()
{
  std::vector<DiscreteFractureConfig> dfms;
  for ( const auto it : m_marker_config )
  {
    const auto & conf = config[it.second];
    DiscreteFractureConfig dfm;
    dfm.region = conf.region;
    dfm.label = it.first;
    dfm.conductivity = conf.conductivity;
    dfm.aperture = conf.aperture;
    dfms.push_back(std::move(dfm));
  }

  return dfms;
}

int EmbeddedFractureManager::find_maximum_face_marker_() const
{
  assert ( m_grid.n_faces() > 0 );
  int max_face_index = std::numeric_limits<int>::lowest();
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
    max_face_index = std::max(max_face_index, face->marker());
  return max_face_index;
}

bool EmbeddedFractureManager::is_fracture(const int face_marker) const
{
  if (m_marker_config.find(face_marker) != m_marker_config.end()) return true;
  else return false;
}

void EmbeddedFractureManager::distribute_mechanical_properties()
{
  find_edfm_cells_and_faces_();

  auto & sda = m_data.sda_data;
  assert( sda.size() == config.size() );
  for (std::size_t i=0; i<config.size(); ++i)
  {
    const size_t n_frac_cells = sda[i].cells.size();
    sda[i].points.assign( n_frac_cells, config[i].body->center() );
    sda[i].dip   .assign( n_frac_cells, config[i].body->plane().dip_angle() );
    sda[i].strike.assign( n_frac_cells, config[i].body->plane().strike_angle() );
    sda[i].cohesion       = config[i].cohesion;
    sda[i].conductivity   = config[i].conductivity;
    sda[i].friction_angle = config[i].friction_angle;
    sda[i].dilation_angle = config[i].dilation_angle;
    sda[i].region         = config[i].region;
  }
}

void EmbeddedFractureManager::
map_mechanics_to_control_volumes(const discretization::DoFNumbering & dofs)
{
  if (m_data.gmcell_to_SDA_flowcells.size() != m_data.geomechanics_grid.n_active_cells())
    m_data.gmcell_to_SDA_flowcells.resize(m_data.geomechanics_grid.n_active_cells());

  for (std::size_t ifrac=0; ifrac<m_data.sda_data.size(); ++ifrac)
  {
    const auto & frac = m_data.sda_data[ifrac];
    for (std::size_t icell=0; icell<frac.cells.size(); ++icell)
    {
      // first map cell cv
      const size_t mech_cell = m_grid.cell(frac.cells[icell]).ultimate_parent().index();
      assert( mech_cell < m_data.gmcell_to_SDA_flowcells.size() );
      for (const size_t iface : frac.faces[icell]){
        if( m_data.coupling[dofs.face_dof(iface)]) // coupled
            m_data.gmcell_to_SDA_flowcells[mech_cell].push_back( dofs.face_dof(iface) );
        else // uncoupled
            config[ifrac].coupled = false;
      }
    }
  }
}

void EmbeddedFractureManager::find_edfm_cells_and_faces_()
{
  std::vector<std::unordered_map<size_t,std::vector<size_t>>> cells_frac_faces(config.size());
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
    if (is_fracture(face->marker()))
    {
      const mesh::Cell *p_neighbor_cell = face->neighbors()[0];
      const mesh::Cell &cell_parent = p_neighbor_cell->ultimate_parent();
      const size_t ifrac = fracture_index_(face->marker());
      cells_frac_faces[ifrac][cell_parent.index()].push_back(face->index());
    }
  auto & sda = m_data.sda_data;
  sda.resize(config.size());
  for (std::size_t ifrac=0; ifrac<sda.size(); ++ifrac)
  {
    const auto & cells_faces = cells_frac_faces[ifrac];
    sda[ifrac].cells.reserve( cells_faces.size() );
    sda[ifrac].faces.reserve( cells_faces.size() );
    for (const auto & it : cells_faces)
    {
      sda[ifrac].cells.push_back(it.first);
      sda[ifrac].faces.push_back(it.second);
    }
  }
}

std::vector<int> EmbeddedFractureManager::get_face_markers() const
{
  std::vector<int> result;
  result.reserve(m_marker_config.size());
  for ( const auto it : m_marker_config )
    result.push_back(it.first);
  return result;
}

// void EmbeddedFractureManager::build_edfm_grid(const discretization::DoFNumbering & dofs)
// {
//   mesh::SurfaceMesh<double> edfm_grid(1e-6);
//   for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
//   {
//     const auto it = m_marker_config.find(face->marker());
//     if (it != m_marker_config.end())
//     {
//       edfm_grid.insert(face->polygon(), it->second);
//       m_data.edfm_cell_mapping.push_back(dofs.face_dof(face->index()));
//     }
//   }

//   m_data.edfm_grid = std::move(edfm_grid);
// }

size_t EmbeddedFractureManager::fracture_index_(const int face_marker) const
{
  return m_marker_config.find(face_marker)->second;
}


}  // end namespace gprs_data
