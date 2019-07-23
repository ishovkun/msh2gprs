#include "OutputDataVTK.hpp"

namespace gprs_data
{

OutputDataVTK::OutputDataVTK(const SimData & sim_data, const mesh::Mesh & grid)
    :
    data(sim_data),
    grid(grid)
{}


void OutputDataVTK::write_output(const std::string & output_path)
{
  save_reservoir_data(output_path + data.config.reservoir_grid_vtk_file);
}

void OutputDataVTK::save_reservoir_mesh(const std::string & fname)
{
  std::ofstream out;
  write_vtk(grid.get_vertices(), grid.get_cells(), grid.shape_ids, out);
  out.close();
}


}
