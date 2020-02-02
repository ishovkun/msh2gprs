#pragma once

#include "PreprocessorConfig.hpp"
#include "SimData.hpp"
#include "discretization/flow/DoFNumbering.hpp"
#include "muparser/muParser.h" // parser for user-defined expressions for reservoir data

namespace gprs_data {

class CellPropertyManager
{
 public:
  CellPropertyManager(const CellPropertyConfig & cell_properties,
                      const std::vector<DomainConfig> & domain_configs,
                      SimData & data);
  // take config and fill out the grid with properties
  void generate_properties();
  // create properties that occured after cell splitting
  void downscale_properties();
  // map coarse geomechanics cells to flow CVs
  void map_mechanics_to_control_volumes(const discretization::DoFNumbering & dofs);
  // remove properties of refined cells after coarsening
  void coarsen_cells();

 private:
  void print_setup_message_();
  void assign_expressions_(const DomainConfig& domain,
                           std::vector<mu::Parser> & parsers,
                           std::vector<double> & vars);
  // return the number of matched cells
  size_t evaluate_expressions_(const DomainConfig& domain,
                               std::vector<mu::Parser> & parsers,
                               std::vector<double> & vars);
  void build_permeability_function_();
  void build_porosity_function_();
  void build_flow_output_property_keys_();

  const CellPropertyConfig        & config;
  const std::vector<DomainConfig> & domains;
  SimData & m_data;
  // number of default variable in config
  // these variables are not output, so I don't save them
  // should be 3=x+y+z
  const size_t m_shift;
  const size_t m_n_unrefined_cells;
};

}  // end namespace gprs_data
