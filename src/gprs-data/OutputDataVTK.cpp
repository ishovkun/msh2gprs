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


void OutputDataVTK::save_reservoir_data(const std::string & fname)
{
  std::ofstream out;
  out.open(fname.c_str());
  IO::VTKWriter::write_geometry(grid.get_vertices(), grid.cells, grid.shape_ids, out);
  IO::VTKWriter::enter_section_cell_data(grid.n_cells(), out);

  // save keywords
  for (std::size_t ivar=0; ivar<data.rockPropNames.size(); ++ivar)
  {
    // if ( data.config.expression_type[ivar] != 1 )  // only mechanics kwds
    //   continue;

    vector<double> property(grid.n_cells());
    const string keyword = data.rockPropNames[ivar];
    for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
      property[cell.index()] = data.vsCellRockProps[cell.index()].v_props[ivar];

    IO::VTKWriter::add_cell_data(property, keyword, out);
  }

  // save multiscale

  out.close();
}


}
