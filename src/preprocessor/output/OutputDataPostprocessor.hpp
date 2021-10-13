#pragma once
#include "SimData.hpp"                          // provides SimData
#include "PreprocessorConfig.hpp"               // provides PreprocessorConfig
#include "discretization/flow/DoFNumbering.hpp" // provides discretization::DoFNumbering
#include "yaml-cpp/yaml.h"                      // provides YAML::Node


namespace gprs_data
{

/**
 * Implements .yaml output for Postrprocessor.py script.
 * Saves a variety of data:
 * 1. Preprocessor output directory.
 * 2. Names of vtk grid files.
 * 3. Mappings matrix/edfm/dfm -> flow dofs
 */
class OutputDataPostprocessor
{
 public:
  /**
   * Constructor
   * 
   * @param  data       : data container
   * @param  config     : container with output file names
   * @param  dofs       : flow dofs
   * @param  output_dir : output path of the preprocessor
   */
  OutputDataPostprocessor(const SimData & data, const PreprocessorConfig & config ,
                          const discretization::DoFNumbering & dofs,
                          const std::string output_dir);
  /**
   * Save config .yaml output file for postprocessor
   * @param  file_name : name of the output file
   */
  void write_output(const std::string file_name);

 private:
  std::vector<size_t> map_matrix_cells_to_control_volumes_();
  const SimData & m_data;
  const PreprocessorConfig & m_config;
  const discretization::DoFNumbering & m_dofs;
  const std::string m_output_dir;
  YAML::Node m_root;
};

}  // end namespace gprs_data
