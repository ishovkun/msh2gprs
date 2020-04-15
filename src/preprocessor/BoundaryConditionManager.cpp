#include "BoundaryConditionManager.hpp"

namespace gprs_data {

BoundaryConditionManager::BoundaryConditionManager(const std::vector<BCConfig> & face_config,
                                                   const std::vector<BCNodeConfig> & node_config,
                                                   SimData & data)
    : _node_config(node_config) , _face_config(face_config),
    _data(data)
{
  build_boundary_conditions_();
}

void BoundaryConditionManager::build_boundary_conditions_()
{
  const auto & grid = _data.geomechanics_grid;
  _node_to_config.resize(grid.n_vertices());
  size_t iface = 0;
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
    if (face->marker() > 0)
    {
      for (size_t iconf=0; iconf<_face_config.size(); ++iconf)
      {
        const auto & conf = _face_config[iconf];
        if ( face->marker() == conf.label )
        {
          if ( conf.type == BoundaryConditionType::dirichlet )
            process_dirichlet_face_(*face, iconf);
          else // if ( conf.type == BoundaryConditionType::neumann )
            process_neumann_face_(conf, iface);
        }
      }
      iface++;
  }
  create_dirichlet_data_();
}

void BoundaryConditionManager::create_dirichlet_data_()
{
  for (size_t v=0; v<_node_to_config.size(); ++v)
    if ( !_node_to_config[v].empty() )
    {
      for (const size_t iconf : _node_to_config[v])
        for (std::size_t i=0; i<3; ++i)
        {
          if ( _face_config[iconf].value[i]  != BCConfig::nan)
          {
            _data.dirichlet_indices[i].push_back(v);
            _data.dirichlet_values[i].push_back( _face_config[iconf].value[i] );
          }
        }
    }
}

void BoundaryConditionManager::process_dirichlet_face_(const mesh::Face & face, const size_t config_index)
{
  for (const size_t v : face.vertices())
    if( std::find( _node_to_config[v].begin(), _node_to_config[v].end(), config_index ) ==
        _node_to_config[v].end())
      _node_to_config[v].push_back( config_index );
}

void BoundaryConditionManager::process_neumann_face_(const BCConfig & conf, const size_t face_index)
{
  _data.neumann_face_indices.push_back( face_index );
  _data.neumann_face_traction.push_back( conf.value );
}

}  // end namespace gprs_data
