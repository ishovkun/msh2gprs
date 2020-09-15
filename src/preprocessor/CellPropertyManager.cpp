#include "CellPropertyManager.hpp"
#include "discretization/flow/DoFNumbering.hpp"
#include "property_expressions/value_functions.hpp"
#include "logger/Logger.hpp"

namespace gprs_data {

CellPropertyManager::
CellPropertyManager(const CellPropertyConfig & cell_properties,
                    const std::vector<DomainConfig> & domain_configs,
                    SimData & data)
    : config(cell_properties), m_data(data), domains(domain_configs),
      m_shift(config.n_default_vars()), m_n_unrefined_cells(data.grid.n_cells_total())
{}

void CellPropertyManager::generate_properties()
{
  print_setup_message_();

  // current number of variables
  const std::size_t n_variables = config.all_vars.size();
  // set up container
  m_data.cell_properties.resize(n_variables);
  for (auto & property : m_data.cell_properties)
    property.assign(m_data.grid.n_cells_total(), 0);

  // save variables name for output
  m_data.property_names.resize(n_variables - m_shift);
  for (std::size_t i=m_shift; i<config.all_vars.size(); ++i)
    m_data.property_names[i - m_shift] = config.all_vars[i];

  // container for evaluated expressions: the properties of a current cell
  std::vector<double> vars(n_variables);

  // loop various domain configs:
  // they may have different number of variables and expressions
  size_t n_matched_cells = 0;
  m_data.coupling.resize( m_data.grid.n_cells_total(), false );
  for (const auto & domain: domains)
  {
    const std::size_t n_expressions = domain.expressions.size();
    // set up muparser
    std::vector<mu::Parser> parsers(n_expressions);
    assign_variables_(parsers, vars);
    assign_expressions_(domain, parsers);
    // run muparser
    n_matched_cells += evaluate_expressions_(domain, parsers, vars);
  }
  if (n_matched_cells != m_data.grid.n_active_cells())
  {
    const std::string msg = "Could not assign properties to all cells. "
                           "Probably missing domain label.";
    throw std::runtime_error(msg);
  }

  std::cout << "Searching for Porosity and Permeability keys" << std::endl;
  build_permeability_function_();
  build_porosity_function_();
  build_volume_mult_function_();
  build_flow_output_property_keys_();
}

size_t CellPropertyManager::evaluate_expressions_(const DomainConfig& domain,
                                                  std::vector<mu::Parser> & parsers,
                                                  std::vector<double> & vars)
{
  const std::size_t n_expressions = domain.expressions.size();
  const std::size_t n_variables = config.all_vars.size();
  const auto & grid = m_data.grid;
  size_t count = 0;  // count number of matched cells
  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
  {
    if (cell->marker() == domain.label) // cells
    {
      count++;
      m_data.coupling[cell->index()] = domain.coupled;
      std::fill(vars.begin(), vars.end(), 0);
      const angem::Point center = cell->center();
      vars[0] = center[0]; // x
      vars[1] = center[1]; // y
      vars[2] = center[2]; // z

      // Evaluate expression -> write into variable
      for (std::size_t i = 0; i < n_expressions; ++i) {
        try {
          vars[domain.local_to_global_vars.at(i)] = parsers[i].Eval();
        } catch (mu::Parser::exception_type &e) {
          const std::string error_msg = "Expression error: " + std::string(e.GetMsg()) +
                                        "\nwhen evaluating expression '" + domain.expressions[i] + "'";
          throw std::runtime_error(error_msg);
        }
      }

      // copy vars to cell properties
      // start from shift to skip x,y,z
      for (std::size_t i = m_shift; i < n_variables; ++i) {
        const size_t property_index = i - m_shift;
        try {
          m_data.cell_properties[property_index][cell->index()] = vars[i];
        } catch (std::out_of_range &e) {  // if subdomain doesn't have property assigned
          throw std::runtime_error( "property not provided for domain");
        }
      }
    } // end match label
  }   // end cell loop
  return count;
}

void CellPropertyManager::assign_variables_(std::vector<mu::Parser> & parsers,
                                            std::vector<double> & vars)
{
  const size_t n_variables = config.all_vars.size();
  for (auto & parser : parsers)
    for (std::size_t j = 0; j < n_variables; ++j)
    {
      try {
        parser.DefineVar(config.all_vars[j], &vars[j]);
      }
      catch (mu::Parser::exception_type &e) {
        const std::string error_msg = "Expression error: " + std::string(e.GetMsg()) +
                                      "\nwhen setting variable '" + config.all_vars[j] + "'";
        throw std::runtime_error(error_msg);
      }
    }
}

void CellPropertyManager::assign_expressions_(const DomainConfig& domain,
                                              std::vector<mu::Parser> & parsers)
{
  const std::size_t n_expressions = domain.expressions.size();
  assign_custom_functions_(parsers);
  for (std::size_t i = 0; i < n_expressions; ++i)
  {
    try {
      parsers[i].SetExpr(domain.expressions[i]);
    } catch (mu::Parser::exception_type &e) {
        const std::string error_msg =
            "Expression error: " + std::string(e.GetMsg()) +
            "\nwhen defining expression '" + domain.expressions[i] + "'";
        throw std::runtime_error(error_msg);
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

void CellPropertyManager::print_setup_message_()
{
  logging::log() << "Setting up function parsers with the following"
                 << " set of variables: "<< std::endl;
  for (std::size_t i=0; i<config.all_vars.size(); ++i)
  {
    logging::log() << config.all_vars[i] << "\t";
    if ((i + 1) % 3 == 0)
      logging::log() << std::endl;
  }
  logging::log() << std::endl;
}

void CellPropertyManager::build_permeability_function_()
{
  bool found_perm_x = false;
  bool found_perm_y = false;
  bool found_perm_z = false;
  for (std::size_t i = 0; i < m_data.property_names.size(); i++)
  {
    const auto & key = m_data.property_names[i];
    if (key == "PERMX")
    {
      found_perm_x = true;
      m_data.permeability_keys[0] = i;
    }
    else if (key == "PERMY")
    {
      found_perm_y = true;
      m_data.permeability_keys[1] = i;
    }
    else if (key == "PERMZ")
    {
      found_perm_z = true;
      m_data.permeability_keys[2] = i;
    }
    else if (key == "PERM")
    {
      // voigh notation
      m_data.permeability_keys[0] = i;
      m_data.permeability_keys[1] = i;
      m_data.permeability_keys[2] = i;
      found_perm_x = true;
      found_perm_y = true;
      found_perm_z = true;
      break;
    }
  }
  if (found_perm_x && found_perm_y && found_perm_z)
        return;
  else
    throw std::invalid_argument("permebility is undefined");
}

void CellPropertyManager::build_porosity_function_()
{
  for (std::size_t i = 0; i < m_data.property_names.size(); i++)
  {
    const auto & key = m_data.property_names[i];
    if (key == "PORO")
    {
      m_data.porosity_key_index = i;
      return;
    }
  }
  throw std::runtime_error("Porosity not provided (Keyword POP)");
}

void CellPropertyManager::build_volume_mult_function_()
{
  for (std::size_t i = 0; i < m_data.property_names.size(); i++)
  {
    const auto & key = m_data.property_names[i];
    if (key == "VFACTOR")
    {
      m_data.vol_mult_index = i;
      return;
    }
  }
}

void CellPropertyManager::build_flow_output_property_keys_()
{
  for (size_t j = 0; j < m_data.property_names.size(); j++)
  {
    if (config.expression_type[j] == ExpressionDomainType::flow)
      m_data.output_flow_properties.push_back(j);
    else if (config.expression_type[j] == ExpressionDomainType::mechanics)
    {
      m_data.output_mech_properties.push_back(j);
      m_data.has_mechanics = true;
    }
  }
}

void CellPropertyManager::map_mechanics_to_control_volumes(const discretization::DoFNumbering & dofs,
                                                           const mesh::Mesh & grid)
{
  m_data.gmcell_to_flowcells.clear();
  m_data.gmcell_to_flowcells.resize(grid.n_active_cells());
  // simdata vector coupled
  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
  {
    for (const auto & conf : domains)
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
  // m_n_unrefined_cells
  for (std::size_t i = 0; i < m_data.property_names.size(); i++)
    m_data.cell_properties[i].erase(m_data.cell_properties[i].begin() + m_n_unrefined_cells,
                                    m_data.cell_properties[i].end());
}

}  // end namespace gprs_data
