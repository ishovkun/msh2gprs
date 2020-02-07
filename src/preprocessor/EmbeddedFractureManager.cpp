#include "EmbeddedFractureManager.hpp"
#include "angem/CollisionGJK.hpp"  // collisionGJK
#include "angem/Collisions.hpp"    // angem::split
#include <utility>                 // provides std::pair
#include <fstream>                 // debug
#include "VTKWriter.hpp"


namespace gprs_data {

using std::vector;
using Point = angem::Point<3,double>;
using std::pair;

EmbeddedFractureManager::
EmbeddedFractureManager(std::vector<EmbeddedFractureConfig> &config,
                        const EDFMMethod edfm_method,
                        SimData & data)
    : config(config), m_method(edfm_method), m_data(data), m_grid(data.grid)
{}

void EmbeddedFractureManager::split_cells()
{
  int face_marker = find_maximum_face_marker_() + 1;
  for (std::size_t i=0; i<config.size(); ++i)
  {
    auto & frac = config[i];
    std::vector<size_t> cells_to_split;
    // iteratively shift fracture if it collides with any grid vertices
    size_t iter = 0;
    while (!find_edfm_cells_(*frac.body, cells_to_split))
    {
      cells_to_split.clear();
      if (++iter > 100)
        throw std::runtime_error("Cannot move fracture to avoid collision with vertices");
    }
    split_cells_(*frac.body, cells_to_split, face_marker);
    m_marker_config.insert({ face_marker, i });
    face_marker++;
  }
  // std::ofstream out;
  // out.open("stuff.vtk");
  // IO::VTKWriter::write_geometry(m_grid, out);
  // out.close();
  // for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
  //   if (face->neighbors().size() > 1)
  // {
  //   std::cout << std::endl;
  //   std::cout << "face->index() = " << face->index() << std::endl;
  //   std::cout << "face->marker() = " << face->marker() << std::endl;
  //   std::cout << "neighbors ";
  //   for (auto pc : face->neighbors())
  //     std::cout << pc->index() << " ";
  //   std::cout << " (";
  //   for (auto pc : face->neighbors())
  //     std::cout << pc->ultimate_parent().index() << " ";
  //   std::cout << ")"<< std::endl;
  //   std::cout << std::endl;
  // }
  // exit(0);
}

void EmbeddedFractureManager::split_cells_(angem::Polygon<double> & fracture,
                                           std::vector<size_t> & cells,
                                           const int face_marker)
{
  const auto & plane = fracture.plane();
  for (const size_t icell : cells)
  {
    mesh::Cell & old_cell = m_grid.cell(icell);
    m_grid.split_cell(old_cell, plane, face_marker);
  }
}

bool EmbeddedFractureManager::find_edfm_cells_(angem::Polygon<double> & fracture,
                                               std::vector<size_t> & cells)
{
  // performs fast collision check
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
        const double dist_vert_center = (polyhedron->center() - vertex).norm();
        if ( std::fabs( fracture.plane().signed_distance(vertex) / dist_vert_center) < 1e-4 )
        {
          const double h = vertices[1].distance(vertices[0]);
          const Point shift = h/5 * fracture.plane().normal();
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
  int max_face_index = m_grid.begin_active_faces()->marker();
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
    sda[i].friction_angle = config[i].friction_angle;
    sda[i].dilation_angle = config[i].dilation_angle;
  }
}

void EmbeddedFractureManager::
map_mechanics_to_control_volumes(const discretization::DoFNumbering & dofs)
{
  if (m_data.gmcell_to_flowcells.size() != m_data.geomechanics_grid.n_active_cells())
    m_data.gmcell_to_flowcells.resize(m_data.geomechanics_grid.n_active_cells());

  for (std::size_t ifrac=0; ifrac<m_data.sda_data.size(); ++ifrac)
  {
    const auto & frac = m_data.sda_data[ifrac];
    for (std::size_t icell=0; icell<frac.cells.size(); ++icell)

    {
      // first map cell cv
      const size_t mech_cell = m_grid.cell(frac.cells[icell]).ultimate_parent().index();
      assert( mech_cell < m_data.gmcell_to_flowcells.size() );
      for (const size_t iface : frac.faces[icell])
        m_data.gmcell_to_flowcells[mech_cell].push_back( dofs.face_dof(iface) );
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

void EmbeddedFractureManager::build_edfm_grid(const discretization::DoFNumbering & dofs)
{
  mesh::SurfaceMesh<double> edfm_grid(1e-6);
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
    if (is_fracture(face->marker()))
    {
      edfm_grid.insert(face->polygon());
      m_data.edfm_cell_mapping.push_back( dofs.face_dof(face->index()) );
    }

  m_data.edfm_grid = std::move(edfm_grid);
}

size_t EmbeddedFractureManager::fracture_index_(const int face_marker) const
{
  return m_marker_config.find(face_marker)->second;
}


}  // end namespace gprs_data
