#include "OutputDataPostprocessor.hpp"
#include <fstream>

namespace gprs_data {

OutputDataPostprocessor::OutputDataPostprocessor(const SimData & data, const PreprocessorConfig & config ,
                                                 const discretization::DoFNumbering & dofs,
                                                 const std::string output_dir)
    : m_data(data), m_config(config), m_dofs(dofs), m_output_dir(output_dir)
{}

void OutputDataPostprocessor::write_output(const std::string file_name)
{
  m_root["output_directory"] = m_config.postprocessor_output_dir;
  m_root["mech_reservoir_grid_file"] = m_output_dir + "/" + m_config.vtk_config.mechanics_reservoir_grid_file;
  m_root["flow_reservoir_grid_file"] = m_output_dir + "/" + m_config.vtk_config.flow_reservoir_grid_file;
  // m_root["edfm_grid_file"] = m_output_dir + "/" + m_config.vtk_config.edfm_grid_file;
  m_root["fracture_grid_file"]  = m_output_dir + "/" + m_config.vtk_config.fracture_grid_file;
  //  map flow vtk cells to flow CVs
  m_root["matrix_cell_to_flow_dof"] = map_matrix_cells_to_control_volumes_();
  m_root["edfm_cell_to_flow_dof"] = m_data.edfm_cell_mapping;
  m_root["dfm_cell_to_flow_dof"] = m_data.dfm_cell_mapping;

  // save file
  std::ofstream out;
  out.open(m_config.postprocessor_file.c_str());
  out << m_root;
  out.close();
}

std::vector<size_t> OutputDataPostprocessor::map_matrix_cells_to_control_volumes_()
{
  const auto & grid = m_data.grid;
  std::vector<size_t> result(grid.n_active_cells());
  size_t cnt = 0;
  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
    result[cnt++] = m_dofs.cell_dof(cell->index());
  return result;
}

}  // end namespace gprs_data
