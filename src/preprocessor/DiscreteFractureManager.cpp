#include "DiscreteFractureManager.hpp"

namespace gprs_data {

DiscreteFractureManager::DiscreteFractureManager(const std::vector<DiscreteFractureConfig> & config,
                                                 SimData & data)
    : m_config(config), m_grid(data.grid), m_data(data)
{}

void DiscreteFractureManager::build_dfm_markers_set_()
{
  for (const auto & frac : m_config)
    m_dfm_markers.insert(frac.label);
}

bool DiscreteFractureManager::is_fracture(const int face_marker) const
{
  if (m_dfm_markers.find(face_marker) != m_dfm_markers.end()) return true;
  else return false;
}

void DiscreteFractureManager::distribute_properties()
{
  std:size_t cv_index;
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
    if (face->neighbors().size() == 2)
      if (is_fracture(face->marker()))
      {
        // find corresponding config
        size_t fracture_index;
        for (size_t i = 0; i < m_config.size(); i++)
          if ( face->marker() == m_config[i].label)
          {
            fracture_index = i;
            break;
          }
        const auto & config = m_config[fracture_index];

        // create a face and fill out the props
        DiscreteFractureFace f;
        f.marker = face->marker();
        f.cv_index = ++cv_index;
        f.coupled = config.coupled;
        f.aperture = config.aperture;
        f.conductivity = config.conductivity;
      }
}

void DiscreteFractureManager::split_faces()
{
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
    if (is_fracture(face->marker()))
      m_grid.mark_for_split(face->index());

  const size_t n_faces_old = m_grid.n_faces();
  const size_t n_vertices_old = m_grid.n_vertices();
  m_data.dfm_grid = m_grid.split_faces();

  if (m_grid.n_vertices() != n_faces_old)
    std::cout << "Split " << m_grid.n_vertices() - n_vertices_old
              << " vertices for dfm fractures"
              << std::endl;
  if (m_grid.n_faces() != n_faces_old)
  {
    std::cout << "Split " << m_grid.n_faces() - n_faces_old
              << " for DFM fractures." << std::endl;
    std::cout << "There was " << n_faces_old << " faces before and "
              << "now there is " << m_grid.n_faces() << " faces."
              << std::endl;
  }
  else
    std::cout << "No DFM faces have been split" << std::endl;
}

}  // end namespace gprs_data
