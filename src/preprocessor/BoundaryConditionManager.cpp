#include "BoundaryConditionManager.hpp"

namespace gprs_data {

BoundaryConditionManager::BoundaryConditionManager(const std::vector<BCConfig> & face_config,
                                                   const std::vector<BCConfig> & node_config,
                                                   SimData & data)
    : _node_config(node_config) , _face_config(face_config), _data(data)
{
  build_boundary_conditions_();
  find_faces_from_expressions_();
  create_dirichlet_data_();
}

void BoundaryConditionManager::build_boundary_conditions_()
{
  const auto & grid = _data.geomechanics_grid;
  _node_to_config.resize(grid.n_vertices());
  size_t iface = 0;
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face, ++iface)
    if (face->marker() > 0)
    {
      for (size_t iconf=0; iconf<_face_config.size(); ++iconf)
      {
        const auto & conf = _face_config[iconf];
        if ( face->marker() == conf.label )
        {
          if (conf.type == BoundaryConditionType::dirichlet)
            process_dirichlet_face_(*face, iconf);
          else // if ( conf.type == BoundaryConditionType::neumann )
            process_neumann_face_(conf, iface);
        }
      }
    }
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

void BoundaryConditionManager::find_faces_from_expressions_()
{
  std::vector<size_t> configs_with_expressions;
  for (size_t i=0; i<_face_config.size(); ++i)
    if (!_face_config[i].expression.empty())
      configs_with_expressions.push_back(i);

  if (configs_with_expressions.empty()) return;

  std::vector<mu::Parser> parsers = create_parsers_(configs_with_expressions);

  const auto & grid = _data.geomechanics_grid;
  size_t iface = 0;
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face, ++iface)
  {
    const auto c = face->center();
    _variables[0] = c[0];
    _variables[1] = c[1];
    _variables[2] = c[2];
    for (size_t i=0; i < configs_with_expressions.size(); ++i)
      if (parsers[i].Eval())
      {
        const auto & conf = _face_config[configs_with_expressions[i]];
        if (conf.type == BoundaryConditionType::neumann)
          process_neumann_face_(conf, iface);
        else
          process_dirichlet_face_(*face, configs_with_expressions[i]);
      }
  }
}

double near(const double x, const double x_target, const double tol)
{
  return std::fabs(x - x_target) < tol;
};

std::vector<mu::Parser> BoundaryConditionManager::create_parsers_(const std::vector<size_t> & configs)
{
  std::vector<mu::Parser> parsers(configs.size());

  for (size_t i=0; i<configs.size(); ++i)
  {
    try {
      parsers[i].DefineVar("X", &_variables[0]);
      parsers[i].DefineVar("Y", &_variables[1]);
      parsers[i].DefineVar("Z", &_variables[2]);
      parsers[i].SetExpr(_face_config[configs[i]].expression);
      parsers[i].DefineFun("near", near);
    } catch (mu::Parser::exception_type &e) {
      const std::string error_msg = "Expression error: " + std::string(e.GetMsg()) +
                                "\nwhen setting variable for boundary conditions";
      throw std::runtime_error(error_msg);
    }
  }
    return parsers;
}

}  // end namespace gprs_data
