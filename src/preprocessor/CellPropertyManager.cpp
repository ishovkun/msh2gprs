#include "CellPropertyManager.hpp"
#include "discretization/flow/DoFNumbering.hpp"
#include "property_expressions/value_functions.hpp"
#include "logger/Logger.hpp"
#include <sstream>

namespace gprs_data {

using namespace std;

CellPropertyManager::
CellPropertyManager(const CellPropertyConfig & cell_properties, SimData & data)
    : _config(cell_properties)
    , m_data(data)
    , m_n_unrefined_cells(data.grid.n_cells_total())
    , _vars(map_variables_())
    , _file_data(read_files_())
{
}

void CellPropertyManager::
generate_properties(std::vector<std::vector<double>> &  properties)
{
  print_setup_message_();
  properties.assign(_vars.size(), vector<double>(m_data.grid.n_cells_total(), 0.f));

  // container for evaluated expressions: the properties of a current cell
  std::vector<double> cell_vars(_vars.size());

  // loop various domain configs:
  // they may have different number of variables and expressions
  size_t n_matched_cells = 0;
  m_data.coupling.resize( m_data.grid.n_cells_total(), false );
  for (const auto & domain: _config.domains)
  {
    const size_t n_expressions = domain.expressions.size();
    // set up muparser
    std::vector<mu::Parser> parsers(n_expressions);
    assign_variables_(parsers, cell_vars);
    assign_expressions_(domain, parsers);
    // run muparser
    n_matched_cells += evaluate_expressions_(domain, parsers, cell_vars, properties);
  }
  if (n_matched_cells != m_data.grid.n_active_cells())
  {
    const std::string msg = "Could not assign properties to all cells. "
                           "Probably missing domain label.";
    throw std::runtime_error(msg);
  }
}

std::unordered_map<std::string, size_t> CellPropertyManager::map_variables_()
{
  size_t nvar = 0;
  std::unordered_map<std::string, size_t> ans;

  for (auto const & var : _config.extra_variables)
    ans[var] = nvar++;

  for (auto const & var : _config.files.variables)
    ans[var] = nvar++;

  for (auto const & domain : _config.domains)
    for (auto const & var : domain.variables)
      ans[var] = nvar++;
  return ans;
}

void CellPropertyManager::evaluate_non_expression_properties_(mesh::Cell const & cell,
                                                              std::vector<double>&cell_vars)
{
  const angem::Point center = cell.center();
  // assign service variables
  cell_vars[_vars[_config.x_kwd]] = center[0];
  cell_vars[_vars[_config.y_kwd]] = center[1];
  cell_vars[_vars[_config.z_kwd]] = center[2];
  cell_vars[_vars[_config.vmult_kwd]] = 1.f;
  // assign file variable values
  for (size_t i = 0; i < _config.files.variables.size(); ++i)
    cell_vars[_vars[_config.files.variables[i]]] = _file_data[i][cell.index()];
}

size_t CellPropertyManager::evaluate_expressions_(const DomainConfig& domain,
                                                  std::vector<mu::Parser> & parsers,
                                                  std::vector<double> & cell_vars,
                                                  std::vector<std::vector<double>> &  properties)
{
  const std::size_t n_expressions = domain.expressions.size();
  const std::size_t n_variables = _vars.size();
  const auto & grid = m_data.grid;
  // we don't want to save x,y,z, and variables from files
  size_t count = 0;  // count number of matched cells
  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
  {
    if (cell->marker() == domain.label) // cells
    {
      count++;
      m_data.coupling[cell->index()] = domain.coupled;
      std::fill(cell_vars.begin(), cell_vars.end(), 0);
      evaluate_non_expression_properties_(*cell, cell_vars);
      // Evaluate expression -> write into variable
      for (std::size_t i = 0; i < n_expressions; ++i) {
        try {
          cell_vars[_vars[domain.variables[i]]] = parsers[i].Eval();
        } catch (mu::Parser::exception_type &e) {
          const std::string error_msg = "Expression error: " + std::string(e.GetMsg()) +
                                        "\nwhen evaluating expression '" + domain.expressions[i] + "'";
          throw std::runtime_error(error_msg);
        }
      }

      // copy local cell_vars to global vectors
      // start from shift to skip x,y,z
      for (std::size_t i = 0; i < _vars.size(); ++i) {
        try {
          properties[i][cell->index()] = cell_vars[i];
        } catch (std::out_of_range &e) {  // if subdomain doesn't have property assigned
          throw std::runtime_error( "property not provided for domain");
        }
      }
    } // end match label
  }   // end cell loop
  return count;
}

std::vector<std::vector<double>> CellPropertyManager::read_files_()
{
  size_t const n_files = _config.files.variables.size();
  size_t const nc = m_n_unrefined_cells;
  std::vector<std::vector<double>> ans(n_files);
  for (size_t i = 0; i < n_files; ++i) {
    ans[i].assign( nc, 0.f );
    std::ifstream in( _config.files.expressions[i], std::ifstream::in );
    size_t cnt = 0;
    while (cnt < nc) {
      if ( !in.good() ) {
        ostringstream msg;
        msg << "File " << _config.files.expressions[i] << " has only " << cnt << " entries. "
            << nc << " entries needed\n";
        throw std::invalid_argument(msg.str());
      }
      in >> ans[i][cnt++];
    }
  }
  return ans;
}

void CellPropertyManager::assign_variables_(std::vector<mu::Parser> & parsers, std::vector<double> & vars)
{
  for (auto & parser : parsers)
    for (auto const & it : _vars) {
      size_t const idx = it.second;
      string const & name = it.first;
      try {
        parser.DefineVar(name, &vars[idx]);
      }
      catch (mu::Parser::exception_type &e) {
        ostringstream msg;
        msg << "Expression error: " << e.GetMsg()
            << "\nwhen setting variable '" << name << "'";
        throw std::runtime_error(msg.str());
      }
    }
}

void CellPropertyManager::assign_expressions_(DomainConfig const& domain, std::vector<mu::Parser> & parsers)
{
  const std::size_t n_expressions = domain.expressions.size();
  assign_custom_functions_(parsers);
  for (size_t i = 0; i < n_expressions; ++i)
  {
    try {
      parsers[i].SetExpr(domain.expressions[i]);
    }
    catch (mu::Parser::exception_type &e) {
      std::ostringstream error_msg;
      error_msg <<  "Expression error: " << e.GetMsg()
                << "\nwhen defining expression '" << domain.expressions[i] << "'";
      throw std::runtime_error(error_msg.str());
    }
  }
}

void CellPropertyManager::assign_custom_functions_(std::vector<mu::Parser> & parsers)
{
  for (auto & parser : parsers)
  {
    parser.DefineFun("clip", property_functions::clip);
  }
}

void CellPropertyManager::print_setup_message_() const
{
  logging::log() << "Setting up function parsers with the following"
                 << " set of variables: "<< std::endl;
  size_t cnt = 0;
  for (auto const & it : _vars)
  {
    logging::log() << it.first << "\t";
    if ((cnt++ + 1) % 5 == 0)
      logging::log() << "\n";
  }
  logging::log() << std::endl;
}

std::vector<size_t> CellPropertyManager::get_custom_flow_keys() const
{
  auto const types = get_property_types();
  std::vector<size_t> ans;
  for (size_t i = 0; i < types.size(); ++i)
    if ( types[i] == VariableType::flow)
      ans.push_back(i);
  return ans;
}
std::vector<size_t> CellPropertyManager::get_custom_mech_keys() const
{
  auto const types = get_property_types();
  std::vector<size_t> ans;
  for (size_t i = 0; i < types.size(); ++i)
    if ( types[i] == VariableType::mechanics)
      ans.push_back(i);
  return ans;
}

std::vector<int> CellPropertyManager::get_permeability_keys() const
{
  vector<int> ans{-1, -1, -1};
  if (_vars.count(_config.perm_kwd))
    std::fill( ans.begin(), ans.end(), _vars[_config.perm_kwd] );
  if (_vars.count(_config.permx_kwd))
    ans[0] = _vars[_config.permx_kwd];
  if (_vars.count(_config.permy_kwd))
    ans[1] = _vars[_config.permy_kwd];
  if (_vars.count(_config.permz_kwd))
    ans[2] = _vars[_config.permz_kwd];

  if ( any_of(ans.begin(), ans.end(), [](auto idx) {return idx < 0;}) ) {
    ostringstream msg;
    msg << "Permeability is undefined. Possible keys: "
        << _config.perm_kwd << " "
        << _config.permx_kwd << " "
        << _config.permy_kwd << " "
        << _config.permz_kwd << std::endl;
    throw std::invalid_argument(msg.str());
  }
  return ans;
}

size_t CellPropertyManager::get_porosity_key() const
{
  if (!_vars.count(_config.poro_kwd)) {
    ostringstream msg;
    msg << "Porosity is not defined. The keyword is : " << _config.poro_kwd;
    throw std::invalid_argument(msg.str());
  }
  return _vars[_config.poro_kwd];
}

void CellPropertyManager::map_mechanics_to_control_volumes(const discretization::DoFNumbering & dofs,
                                                           const mesh::Mesh & grid)
{
  m_data.gmcell_to_flowcells.clear();
  m_data.gmcell_to_flowcells.resize(grid.n_active_cells());
  // simdata vector coupled
  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
  {
    for (const auto & conf : _config.domains)
      if (conf.label == cell->marker())
        if (conf.coupled)
        {
          m_data.gmcell_to_flowcells[m_data.mech_numbering->cell_dof(cell->index())]
              .push_back(dofs.cell_dof(cell->index()));
        }
  }
}

void CellPropertyManager::downscale_properties()
{
  const auto & grid = m_data.grid;
  const size_t n_raw_cells = grid.n_cells_total();
  for (std::size_t i = 0; i < m_data.property_names.size(); i++)
    m_data.cell_properties[i].resize( n_raw_cells );

  m_data.coupling.resize( n_raw_cells );
  for (auto raw = grid.begin_cells() + m_n_unrefined_cells; raw != grid.end_cells(); ++raw )
  {
    for (std::size_t i = 0; i < m_data.property_names.size(); i++)
      m_data.cell_properties[i][raw->index()] = m_data.cell_properties[i][raw->ultimate_parent().index()];
    m_data.coupling[raw->index()] = m_data.coupling[raw->ultimate_parent().index()];
  }
}

void CellPropertyManager::coarsen_cells()
{
  for (std::size_t i = 0; i < m_data.property_names.size(); i++)
    m_data.cell_properties[i].erase(m_data.cell_properties[i].begin() + m_n_unrefined_cells,
                                    m_data.cell_properties[i].end());
}

std::vector<std::string> CellPropertyManager::get_property_names() const
{
  size_t const n_output_vars = _vars.size();
  vector<string> ans(n_output_vars);
  for (auto const & domain : _config.domains)
    for (auto const & var : domain.variables)
      ans[_vars[var]] = var;
  return ans;
}

std::vector<VariableType> CellPropertyManager::get_property_types() const
{
  size_t const n_output_vars = _vars.size();
  vector<VariableType> ans(n_output_vars);
  for (auto const & domain : _config.domains)
    for (auto const & var : domain.variables)
      ans[_vars[var]] = domain.type;
  for (auto const & var : _config.extra_variables)
      ans[_vars[var]] = VariableType::service;
  for (auto const & var : _config.files.variables)
      ans[_vars[var]] = VariableType::service;
  return ans;
}


}  // end namespace gprs_data
