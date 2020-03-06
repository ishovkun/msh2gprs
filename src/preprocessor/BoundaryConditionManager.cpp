#include "BoundaryConditionManager.hpp"

namespace gprs_data {

BoundaryConditionManager::BoundaryConditionManager(const std::vector<BCConfig> & face_config,
                                                   const std::vector<BCNodeConfig> & node_config,
                                                   SimData & data)
    : m_node_config(node_config) , m_face_config(face_config),
    m_data(data)
{
}

void BoundaryConditionManager::create_properties()
{
  const std::unordered_map<size_t,size_t> face_to_frac = create_boundary_faces_();
  create_neumann_faces_(face_to_frac);
}

std::unordered_map<size_t,size_t> BoundaryConditionManager::create_boundary_faces_() const
{
  std::map<int,size_t> marker_config;
  for (std::size_t i=0; i<m_face_config.size(); ++i)
    marker_config.insert({m_face_config[i].label, i});

  const auto & grid = m_data.grid;
  std::unordered_map<size_t, size_t> face_config;
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
  {
    const auto it = marker_config.find(face->marker());
    if (it != marker_config.end())
      face_config[face->index()] = it->second;
  }
  return face_config;
}

void BoundaryConditionManager::create_neumann_faces_(const std::unordered_map<size_t,size_t> & face_to_frac)
{
  for (const auto & it : face_to_frac)
  {
    if (m_face_config[it.second].type == BoundaryConditionType::neumann)
    {
      m_data.neumann_face_indices.push_back( it.first );
      m_data.neumann_face_traction.push_back(m_face_config[it.second].value);
    }
  }
}

}  // end namespace gprs_data
