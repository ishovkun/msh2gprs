#include "EmbeddedFractureManager.hpp"
#include "angem/CollisionGJK.hpp"  // collisionGJK
#include "angem/Collisions.hpp"    // angem::split
#include <utility>                 // provides std::pair

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
  for (auto & frac : config)  // non-const since we can shift it
  {
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
    m_edfm_markers.insert(face_marker);
    face_marker++;
  }
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
    // const angem::Polyhedron<double> & poly_cell = *pol;
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
  size_t i = 0;
  for (const int marker : m_edfm_markers)
  {
    const auto & conf = config[i];
    DiscreteFractureConfig dfm;
    dfm.label = marker;
    dfm.conductivity = conf.conductivity;
    dfm.aperture = conf.aperture;
    dfms.push_back(std::move(dfm));
    i++;
  }

  return dfms;
}

int EmbeddedFractureManager::find_maximum_face_marker_() const
{
  int max_face_index = m_grid.begin_active_faces()->marker();
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
    max_face_index = std::max(max_face_index, face->marker());
  return max_face_index;
}

bool EmbeddedFractureManager::is_fracture(const int face_marker) const
{
  if (m_edfm_markers.find(face_marker) != m_edfm_markers.end()) return true;
  else return false;
}

void EmbeddedFractureManager::distribute_mechanical_properties()
{
  auto & sda = m_data.sda_data;
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

void EmbeddedFractureManager::map_mechanics_to_control_volumes()
{
  // auto & sda = m_data.sda_data;
  // // first map edfm face to pair frac-index
  // std::unordered_map<size_t, pair<size_t, size_t>> cell_to_frac;
  // for (std::size_t ifrac=0; ifrac<sda.size(); ++ifrac)
  //   for (std::size_t icell=0; icell<sda[ifrac].cells.size(); ++icell)
  //     cell_to_frac[sda[ifrac][icell]] = {ifrac, icell};

  // for (std::size_t i=0; i<config.size(); ++i)
  //   if (!sda[i].cells.empty())
  //     sda[i].cvs.resize(sda[i].cells.size());

  // for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
  //   if (is_fracture(face->marker()))
  //   {
  //     const mesh::Cell *p_neighbor_cell = face->neighbors()[0];
  //     const mesh::Cell &cell_parent = p_neighbor_cell->ultimate_parent();
  //     const auto it_pair_ifrac_icell = cell_to_frac.find(cell_parent.index());
  //     const size_t ifrac = it_pair_ifrac_icell->second.first;
  //     const size_t icell = it_pair_ifrac_icell->second.second;
  //     // sda[ifrac].cvs[icell].push_back();
  //   }

  //   //   }
  //   // }
}

std::vector<int> EmbeddedFractureManager::get_face_markers() const
{
  return std::vector<int>(m_edfm_markers.begin(), m_edfm_markers.end());
}

void EmbeddedFractureManager::build_edfm_grid()
{
  mesh::SurfaceMesh<double> edfm_grid(1e-6);
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
    if (is_fracture(face->marker()))
      edfm_grid.insert(face->polygon());
  m_data.edfm_grid = std::move(edfm_grid);
}

}  // end namespace gprs_data
