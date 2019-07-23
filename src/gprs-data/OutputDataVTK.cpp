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
  if (data.dfm_faces.size() > 0)
    save_dfm_data(output_path + data.config.discrete_frac_file);
  if (!data.vEfrac.empty())
    save_edfm_data(output_path + data.config.efrac_file);
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
  const auto & ms = data.ms_data;
  if (!ms.partitioning.empty())
  {
    IO::VTKWriter::add_cell_data(ms.partitioning, "partitioning", out);

    // support regions
    std::vector<int> support_value(grid.n_cells());
    for (std::size_t block = 0; block < ms.n_blocks; block++)
    {
      for (size_t cell = 0; cell < grid.n_cells(); ++cell)
      {
        if (cell == ms.centroids[block])
          support_value[cell] = 3;
        else if (ms.support_boundary_cells[block].find(cell) != ms.support_boundary_cells[block].end())
          support_value[cell] = 2;
        else if (ms.support_internal_cells[block].find(cell) != ms.support_internal_cells[block].end())
          support_value[cell] = 1;
        else
          support_value[cell] = 0;
      }

      IO::VTKWriter::add_cell_data(support_value, "support-" + std::to_string(block), out);
    }
  }

  out.close();
}


void OutputDataVTK::save_dfm_data(const std::string & fname)
{
  std::cout << "Saving DFM mesh file: " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());
  IO::VTKWriter::write_surface_geometry(data.dfm_master_grid.get_vertices(),
                                        data.dfm_master_grid.get_polygons(),
                                        out);
  out.close();
}


void OutputDataVTK::save_edfm_data(const std::string & fname)
{
  std::cout << "Saving EDFM mesh file: " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());

  // write vtk data
  std::size_t n_efrac_vertices = 0;
  // make up a vector of all sda vertices
  for (const auto & efrac : data.vEfrac)
    n_efrac_vertices += efrac.mesh.n_vertices();

  std::vector<angem::Point<3,double>> efrac_verts(n_efrac_vertices);
  std::size_t ivertex = 0;
  for (const auto & efrac : data.vEfrac)
    for (const auto & p : efrac.mesh.get_vertices())
    {
      efrac_verts[ivertex] = p;
      ivertex++;
    }

  std::size_t n_efrac_elements = 0;
  std::size_t vind_size_total = 0;

  for (const auto & efrac : data.vEfrac)
  {
    n_efrac_elements += efrac.mesh.n_polygons();
    for (const auto & poly : efrac.mesh.get_polygons())
      vind_size_total += poly.size();
  }

  std::vector<std::vector<std::size_t>> efrac_cells(n_efrac_elements);
  std::size_t ielement = 0;
  std::size_t shift = 0;
  for (const auto & efrac : data.vEfrac)
  {
    for (const auto & cell : efrac.mesh.get_polygons())
    {
      efrac_cells[ielement].resize(cell.size());
      for (short v=0; v<cell.size(); ++v)
        efrac_cells[ielement][v] = shift + cell[v];
      ielement++;
    }
    shift += efrac.mesh.n_vertices();
  }

  IO::VTKWriter::write_surface_geometry(efrac_verts, efrac_cells, out);


  out.close();
}


}
