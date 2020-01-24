#include "DiscreteFractureManager.hpp"

namespace gprs_data {

DiscreteFractureManager::DiscreteFractureManager(const std::vector<DiscreteFractureConfig> & config,
                                                 SimData & data)
    : m_config(config), m_grid(data.grid), m_data(data)
{
  build_dfm_markers_set_();
}

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
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
      if (is_fracture(face->marker()))
      {
        if (face->neighbors().size() < 2)
          throw std::runtime_error("degenerate grid. Fracture face " +
                                   std::to_string(face->index()) + " marker " +
                                   std::to_string(face->marker()) + " " +
                                   std::to_string(face->neighbors().size()) +
                                   " neighbors ");

        // find corresponding config
        size_t fracture_index = std::numeric_limits<size_t>::max();
        for (size_t i = 0; i < m_config.size(); i++)
          if ( face->marker() == m_config[i].label)
          {
            fracture_index = i;
            break;
          }
        assert( fracture_index != std::numeric_limits<size_t>::max() );
        const auto & config = m_config[fracture_index];

        // create a face and fill out the props
        DiscreteFractureFace f;
        f.marker = face->marker();
        f.coupled = config.coupled;
        f.aperture = config.aperture;
        f.conductivity = config.conductivity;
        // compute custom data as weighted average of the bounding cells
        const auto cells = face->neighbors();
        const auto cell1 = *cells[0];
        const auto cell2 = *cells[1];
        f.custom_flow_data.resize(m_data.output_flow_properties.size());
        for (size_t j = 0; j < m_data.output_flow_properties.size(); ++j)
        {
          const size_t key = m_data.output_flow_properties[j];
          f.custom_flow_data[j] +=
              (m_data.cell_properties[key][cell1.index()] * cell1.volume() +
               m_data.cell_properties[key][cell2.index()] * cell2.volume() ) /
              ( cell1.volume() + cell2.volume() ) ;
        }
        m_data.dfm_faces.insert({face->index(), std::move(f)});
      }
}

void DiscreteFractureManager::split_faces(mesh::Mesh & grid)
{
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
    if (is_fracture(face->marker()))
      grid.mark_for_split(face->index());

  const size_t n_faces_old = grid.n_faces();
  const size_t n_vertices_old = grid.n_vertices();

  grid.split_faces();

  if (grid.n_vertices() != n_faces_old)
    std::cout << "Split " << grid.n_vertices() - n_vertices_old
              << " vertices for dfm fractures"
              << std::endl;
  if (grid.n_faces() != n_faces_old)
  {
    std::cout << "Split " << grid.n_faces() - n_faces_old
              << " for DFM fractures." << std::endl;
    std::cout << "There was " << n_faces_old << " faces before and "
              << "now there is " << grid.n_faces() << " faces."
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

size_t DiscreteFractureManager::count_dfm_faces() const
{
  size_t n_dfm_faces = 0;
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
    if (face->neighbors().size() == 2)
      if (is_fracture(face->marker()))
        n_dfm_faces++;
  return n_dfm_faces;
}

std::vector<int> DiscreteFractureManager::get_face_markers() const
{
  return std::vector<int>(m_dfm_markers.begin(), m_dfm_markers.end());
}

mesh::SurfaceMesh<double> DiscreteFractureManager::build_dfm_grid(const mesh::Mesh & grid) const
{
  mesh::SurfaceMesh<double> dfm_grid(1e-6);
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
    if (face->neighbors().size() == 2)
      if (is_fracture(face->marker()))
        dfm_grid.insert(face->polygon());

  return dfm_grid;
}

std::vector<size_t> DiscreteFractureManager::map_dfm_grid_to_flow_dofs(const mesh::Mesh & grid,
                                                                       const discretization::DoFNumbering & dofs) const
{
  std::vector<size_t> result;
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
    if (face->neighbors().size() == 2)
      if (is_fracture(face->marker()))
        result.push_back( dofs.face_dof(face->index()) );
  return result;
}


}  // end namespace gprs_data
