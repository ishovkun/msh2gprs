#include "EmbeddedFractureManager.hpp"
#include "angem/CollisionGJK.hpp"  // collisionGJK
#include "angem/Collisions.hpp"    // angem::split

namespace gprs_data {

using std::vector;
using Point = angem::Point<3,double>;

EmbeddedFractureManager::
EmbeddedFractureManager(const std::vector<EmbeddedFractureConfig> &config,
                        const EDFMMethod edfm_method,
                        SimData & data)
    : config(config), m_method(edfm_method), data(data), m_grid(data.grid)
{}

void EmbeddedFractureManager::split_cells()
{
  int face_marker = find_maximum_face_marker_() + 1;
  for (auto & frac : config)  // non-const since we can shift it
  {
    vector<size_t> cells;
    // iteratively shift fracture if it collides with any grid vertices
    size_t iter = 0;
    while (!find_edfm_cells_(*frac.body, cells))
    {
      cells.clear();
      if (++iter > 100)
        throw std::runtime_error("Cannot move fracture to avoid collision with vertices");
    }

    split_cells_(*frac.body, cells, face_marker);
    std::cout << "should be new face marker " << face_marker << std::endl;
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
    // std::cout << "splitting " << old_cell.index() << std::endl;
    m_grid.split_cell(old_cell, plane, face_marker);
  }
}

bool EmbeddedFractureManager::find_edfm_cells_(angem::Polygon<double> & fracture,
                                               std::vector<size_t> & cells) const
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

  std::cout << "max_face_index = " << max_face_index << std::endl;
  return max_face_index;
}

bool EmbeddedFractureManager::is_fracture(const int face_marker) const
{
  if (m_edfm_markers.find(face_marker) != m_edfm_markers.end()) return true;
  else return false;
}

std::vector<discretization::ControlVolumeData> EmbeddedFractureManager::
extract_control_volume_data(const std::vector<discretization::ControlVolumeData> & cv_data,
                            const size_t n_dfm_faces)
{
  std::vector<discretization::ControlVolumeData> edfm_cv_data;
  size_t n_edfm_faces = 0;
  for (std::size_t i=n_dfm_faces; i<cv_data.size(); ++i)
    if ( cv_data[i].type == discretization::ControlVolumeType::face)
      n_edfm_faces++;
  edfm_cv_data.reserve(n_edfm_faces);
  for (std::size_t i=n_dfm_faces; i<n_dfm_faces + n_edfm_faces; ++i)
    if ( cv_data[i].type == discretization::ControlVolumeType::face)
      edfm_cv_data.push_back( cv_data[i] );
  return edfm_cv_data;
}

void EmbeddedFractureManager::
extract_flow_data(const std::vector<discretization::ControlVolumeData> & cv_data,
                  const std::vector<discretization::ConnectionData> & con_data,
                  const size_t n_dfm_faces, const size_t n_cells)
{
  // count edfm entries
  const size_t min_edfm_index = n_dfm_faces;
  size_t n_edfm_faces = 0;
  for (std::size_t i=n_dfm_faces; i<cv_data.size(); ++i)
    if ( cv_data[i].type == discretization::ControlVolumeType::face)
      n_edfm_faces++;
  std::cout << "n_edfm_faces = " << n_edfm_faces << std::endl;

  const size_t max_edfm_index = min_edfm_index + n_edfm_faces - 1;
  std::cout << "here it comes: " << min_edfm_index << " " << max_edfm_index << std::endl;
  for (const auto & con: con_data)
    if (con.type == discretization::ConnectionType::matrix_fracture)
    {
      assert( con.elements.size() == 2 );
      // std::cout << "here " << con.elements[0] << " " << con.elements[1] << std::endl;
      // identify which of connection elements is edfm
      size_t frac_element = con.elements.size();
      if      ( min_edfm_index <= con.elements[0] &&
                max_edfm_index >= con.elements[0])
        frac_element = 0;
      else if ( min_edfm_index <= con.elements[1] &&
                max_edfm_index >= con.elements[1])
        frac_element = 1;
      if (frac_element == 2) continue;

      std::cout << con.elements[frac_element] << std::endl;
      // connect to cell parent since we split edfm cells

    }
    else if (con.type == discretization::ConnectionType::fracture_fracture)
    {
      
    }

  assert ( false && "write extraction code" );
}

}  // end namespace gprs_data
