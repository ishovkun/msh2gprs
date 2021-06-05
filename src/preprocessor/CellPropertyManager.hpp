#pragma once

#include "PreprocessorConfig.hpp"
#include "SimData.hpp"
#include "discretization/flow/DoFNumbering.hpp"
#include "muparser/muParser.h" // parser for user-defined expressions for reservoir data

namespace gprs_data {

/* This class implementes a manager for grid cell properties.
 * Its main responsibiilities are to 
 * (1) distribute properties across the grid and grid subdomains.
 * (2) upscale and downscale properties upon grid refinement.  */
class CellPropertyManager
{
 public:
  /**
   * Constructor. 
   * Requires a generic config and config for each subdomain.
   * 
   * @param  {CellPropertyConfig} cell_properties       : generic config
   * @param  {SimData} data                             : container for output data
   */
  CellPropertyManager(const CellPropertyConfig & cell_properties, SimData & data);
  // take config and fill out the grid with properties
  void generate_properties(std::vector<std::vector<double>> & cell_properties);
  // get names of properties
  std::vector<std::string> get_property_names() const;
  // get types of properties
  std::vector<VariableType> get_property_types() const;
  // create properties that occured after cell splitting
  void downscale_properties();
  // map coarse geomechanics cells to flow CVs
  void map_mechanics_to_control_volumes(const discretization::DoFNumbering & dofs,
                                        const mesh::Mesh & grid);
  // remove properties of refined cells after coarsening
  void coarsen_cells();
  // get indices of permeabilities in properties vector
  std::vector<int> get_permeability_keys() const;
  size_t get_porosity_key() const;
  size_t get_volume_mult_key() const {return _vars[_config.vmult_kwd];}

 private:
  void print_setup_message_() const;
  void assign_variables_(std::vector<mu::Parser> & parsers, std::vector<double> & vars);
  void assign_expressions_(const DomainConfig& domain, std::vector<mu::Parser> & parsers);
  // return the number of matched cells
  size_t evaluate_expressions_(const DomainConfig& domain,
                               std::vector<mu::Parser> & parsers,
                               std::vector<double> & cell_vars,
                               std::vector<std::vector<double>> &  cell_properties);
  void assign_custom_functions_(std::vector<mu::Parser> & parsers);
  void build_permeability_function_();

  void build_flow_output_property_keys_();
  // void build_volume_mult_function_();
  void evaluate_non_expression_properties_(mesh::Cell const & cell, std::vector<double>&cell_vars);
  std::unordered_map<std::string, size_t> map_variables_();
  std::vector<std::vector<double>> read_files_();

  CellPropertyConfig const & _config;
  SimData & m_data;
  const size_t m_n_unrefined_cells;
  mutable std::unordered_map<std::string, size_t> _vars;
  std::vector<std::vector<double>> _file_data;
};

}  // end namespace gprs_data
