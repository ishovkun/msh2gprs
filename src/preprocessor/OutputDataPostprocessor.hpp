#pragma once
#include "SimData.hpp"
#include "PreprocessorConfig.hpp"
#include "discretization/flow/DoFNumbering.hpp"
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep


namespace gprs_data
{

class OutputDataPostprocessor
{
 public:
  OutputDataPostprocessor(const SimData & data, const PreprocessorConfig & config ,
                          const discretization::DoFNumbering & dofs,
                          const std::string output_dir);
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
