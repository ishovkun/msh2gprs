#include "BoundaryConditionManager.hpp"
#include "bc_expressions/location_functions.hpp"
#include "bc_expressions/value_functions.hpp"

namespace gprs_data {

BoundaryConditionManager::BoundaryConditionManager(const std::vector<BCConfig> & face_config,
                                                   const std::vector<BCConfig> & node_config,
                                                   SimData & data)
    : _node_config(node_config) , _face_config(face_config), _data(data)
{
  // create value parsers
  create_value_parsers_(_face_config, _face_value_parsers);
  create_value_parsers_(_node_config, _node_value_parsers);

  build_boundary_conditions_();
  find_faces_from_expressions_();
  create_dirichlet_data_();
}

void BoundaryConditionManager::build_boundary_conditions_()
{
  const auto & grid = _data.geomechanics_grid;
  _node_to_config.resize(grid.n_vertices());
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
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
            process_neumann_face_(*face, iconf);

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
      {
        const auto coord = _data.geomechanics_grid.vertex(v);
        _variables[0] = coord[0];
        _variables[1] = coord[1];
        _variables[2] = coord[2];
        _variables[3] = BCConfig::nan;
        for (std::size_t i=0; i<3; ++i)
        {
          const double nan_value = BCConfig::nan;
          double value;
          try {value = _face_value_parsers[iconf][i].Eval();}
          catch (mu::Parser::exception_type &e) {
            const std::string error_msg = "Expression error: " + std::string(e.GetMsg()) +
                "\nwhen setting evaluating boundary condition";
            throw std::runtime_error(error_msg);
          }
          if (value != nan_value)
          {
            _data.dirichlet_indices[i].push_back(v);
            _data.dirichlet_values[i].push_back(value);
          }
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

void BoundaryConditionManager::process_neumann_face_(const mesh::Face & face, const size_t config_index)
{
  _data.neumann_face_indices.push_back( face.index() );
  const auto c = face.center();
  _variables[0] = c[0];
  _variables[1] = c[1];
  _variables[2] = c[2];
  _variables[3] = BCConfig::nan;
  angem::Point<3,double> value;
  for (size_t i=0; i<3; ++i)
    value[i] = _face_value_parsers[config_index][i].Eval();
  _data.neumann_face_traction.push_back( value );
}

void BoundaryConditionManager::find_faces_from_expressions_()
{
  std::vector<size_t> configs_with_expressions;
  for (size_t i=0; i<_face_config.size(); ++i)
    if (!_face_config[i].location_expression.empty())
      configs_with_expressions.push_back(i);

  if (configs_with_expressions.empty()) return;

  std::vector<mu::Parser> parsers = create_location_parsers_(configs_with_expressions);

  const auto & grid = _data.geomechanics_grid;
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
  {
    const auto c = face->center();
    _variables[0] = c[0];
    _variables[1] = c[1];
    _variables[2] = c[2];
    _variables[3] = BCConfig::nan;
    for (size_t i=0; i < configs_with_expressions.size(); ++i)
      if (parsers[i].Eval())  // found location specified by location function
      {
        const auto & conf = _face_config[configs_with_expressions[i]];
        if (conf.type == BoundaryConditionType::neumann)
          process_neumann_face_(*face, configs_with_expressions[i]);
        else
          process_dirichlet_face_(*face, configs_with_expressions[i]);
      }
  }
}

std::vector<mu::Parser> BoundaryConditionManager::create_location_parsers_(const std::vector<size_t> & configs)
{
  std::vector<mu::Parser> parsers(configs.size());

  for (size_t i=0; i<configs.size(); ++i)
  {
    try {
      parsers[i].DefineVar("X", &_variables[0]);
      parsers[i].DefineVar("Y", &_variables[1]);
      parsers[i].DefineVar("Z", &_variables[2]);
      parsers[i].DefineVar("nan", &_variables[3]);
      parsers[i].SetExpr(_face_config[configs[i]].location_expression);
      parsers[i].DefineFun("almost_equal", bc_functions::almost_equal);
      parsers[i].DefineFun("near", bc_functions::near);
    }
    catch (mu::Parser::exception_type &e) {
      const std::string error_msg = "Expression error: " + std::string(e.GetMsg()) +
                                "\nwhen setting variable for boundary conditions";
      throw std::runtime_error(error_msg);
    }
  }
    return parsers;
}

void BoundaryConditionManager::create_value_parsers_(const std::vector<BCConfig> & config,
                                                     std::vector<std::array<mu::Parser,3>> & parsers)
{
  parsers.resize(config.size());
  for (size_t i=0; i<config.size(); ++i)
    for (size_t j=0; j<3; ++j)
    {
      try {
        parsers[i][j].DefineVar("X", &_variables[0]);
        parsers[i][j].DefineVar("Y", &_variables[1]);
        parsers[i][j].DefineVar("Z", &_variables[2]);
        parsers[i][j].DefineVar("nan", &_variables[3]);
        parsers[i][j].DefineFun("cantilever_beam_end_shear_ux",
                                bc_functions::cantilever_beam_end_shear_ux);
        parsers[i][j].DefineFun("cantilever_beam_end_shear_uy",
                                bc_functions::cantilever_beam_end_shear_uy);
        parsers[i][j].DefineFun("cantilever_beam_end_shear_uz",
                                bc_functions::cantilever_beam_end_shear_uz);
        if (config[i].values_expressions[j].empty())
        {
          throw std::invalid_argument("expression " + std::to_string(i) + " is empty");
        }
        parsers[i][j].SetExpr(config[i].values_expressions[j]);
      }
    catch (mu::Parser::exception_type &e) {
      const std::string error_msg = "Expression error: " + std::string(e.GetMsg()) +
                                "\nwhen setting variable for boundary conditions";
      throw std::runtime_error(error_msg);
    }
  }
}

}  // end namespace gprs_data
