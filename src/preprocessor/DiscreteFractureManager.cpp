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
  m_data.dfm_faces.clear();
  std:size_t cv_index = 0;
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
        m_data.dfm_faces.insert({face->index(), std::move(f)});
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

std::vector<DiscreteFractureConfig>
DiscreteFractureManager::combine_configs(const std::vector<DiscreteFractureConfig> & config1,
                                         const std::vector<DiscreteFractureConfig> & config2)
{
  std::vector<DiscreteFractureConfig> combined;
  combined.reserve( config1.size() + config2.size() );
  for (const auto & conf : config1)
    combined.push_back(conf);
  for (const auto & conf : config2)
    combined.push_back(conf);
  return combined;
}

void DiscreteFractureManager::build_reservoir_cell_numbering()
{
  size_t cv_index = 0;
  // distribute dfm indices
  for (auto & it_face: m_data.dfm_faces)
    it_face.second.cv_index = ++cv_index;
  // distribute cell indices
  m_data.cell_cv_indices.clear();
  m_data.cell_cv_indices.reserve( m_grid.n_cells() );
  for (auto cell = m_grid.begin_active_cells(); cell != m_grid.end_active_cells(); ++cell)
    m_data.cell_cv_indices.push_back( ++cv_index );
}

}  // end namespace gprs_data
