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
    save_dfm_data(output_path + data.config.dfm_grid_vtk_file);
  if (!data.vEfrac.empty())
    save_edfm_data(output_path + data.config.edfm_grid_vtk_file);
  if (!data.well_vertices.points.empty())
    save_well_data(output_path + data.config.wells_vtk_file);
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
    vector<double> property(grid.n_cells());
    const string keyword = data.rockPropNames[ivar];
    for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
      property[cell.index()] = data.vsCellRockProps[cell.index()].v_props[ivar];

    IO::VTKWriter::add_data(property, keyword, out);
  }

  // save multiscale flow data
  if (!data.ms_flow_data.partitioning.empty())
  {
    IO::VTKWriter::add_data(data.ms_flow_data.partitioning, "partitioning-flow", out);
    saveMultiScaleSupport(data.ms_mech_data, grid.n_cells(), "support-flow-", out);
  }
  if (!data.ms_mech_data.partitioning.empty())
  {
    const auto & ms = data.ms_mech_data;
    IO::VTKWriter::add_data(ms.partitioning, "partitioning-mech", out);
    IO::VTKWriter::enter_section_point_data(grid.n_vertices(), out);
    // saveMultiScaleSupport(msm, grid.n_vertices(), "support-mech-", out);

    // I'm just gonna save vertices and boundaries
    std::vector<int> output(grid.n_vertices(), 0);
    for (std::size_t coarse = 0; coarse < ms.n_coarse; coarse++)
    {
      for (const size_t i : ms.support_boundary[coarse])
        output[i] = 1;
    }
    for (const size_t i : ms.centroids)
      output[i] = 2;
    IO::VTKWriter::add_data(output, "support-mech", out);
  }
  out.close();
}


void OutputDataVTK::save_dfm_data(const std::string & fname)
{
  std::cout << "Saving DFM mesh file: " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());
  IO::VTKWriter::write_surface_geometry(data.dfm_master_grid.get_vertices(),
                                        data.dfm_master_grid.get_polygons(), out);
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


void OutputDataVTK::saveMultiScaleSupport(const multiscale::MultiScaleOutputData & ms,
                                          const std::size_t                        size,
                                          const std::string                      & prefix,
                                          std::ofstream                          & out)
{
  // coarse nodes
  std::vector<int> support_value(size);
  for (std::size_t coarse = 0; coarse < ms.n_coarse; coarse++)
  {
    for (std::size_t i=0; i<size; ++i)
    {
      if (i == ms.centroids[coarse])
        support_value[i] = 3;
      else if (ms.support_boundary[coarse].find(i) != ms.support_boundary[coarse].end())
        support_value[i] = 2;
      else if (ms.support_internal[coarse].find(i) != ms.support_internal[coarse].end())
        support_value[i] = 1;
      else
        support_value[i] = 0;
    }

    IO::VTKWriter::add_data(support_value, prefix + std::to_string(coarse), out);
  }  // end coarse loop
}


void OutputDataVTK::save_well_data(const std::string & fname)
{
  // vtk well geometry
  IO::VTKWriter::write_well_trajectory(data.well_vertices.points,
                                       data.well_vertex_indices,
                                       fname);
}

}
