#include "OutputDataVTK.hpp"
#include "logger/Logger.hpp"

namespace gprs_data
{

OutputDataVTK::OutputDataVTK(const SimData & data,
                             const VTKOutputConfig config)
    :
    _data(data),
    m_flow_grid(data.grid),
    m_mech_grid(data.geomechanics_grid),
    m_config(config)
{}


void OutputDataVTK::write_output(const std::string & output_path) const
{
  save_reservoir_flow_data_(output_path + "/" + m_config.flow_reservoir_grid_file);
  save_reservoir_mechanics_data_(output_path + "/" + m_config.mechanics_reservoir_grid_file);
  if (!_data.fracture_grid.empty())
    save_dfm_data(output_path + "/"+ m_config.fracture_grid_file);
  // if (!_data.edfm_grid.empty())
  //   save_edfm_data(output_path + "/" + m_config.edfm_grid_file);
  if (!_data.wells.empty())
    save_wells_(output_path + "/"+ m_config.wells_file);
}


void OutputDataVTK::save_reservoir_flow_data_(const std::string & fname) const
{
  logging::log() << "writing " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());
  IO::VTKWriter::write_geometry(m_flow_grid, out);
  IO::VTKWriter::enter_section_cell_data(m_flow_grid.n_active_cells(), out);

  const auto & grid = m_flow_grid;
  // save grid markers
  {
    std::vector<double> property(grid.n_active_cells());
    const std::string keyword = "Marker";
    size_t cell_index = 0;
    for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
      property[cell_index++] = cell->marker();
    IO::VTKWriter::add_data(property, keyword, out);
  }
  // save porosity
  {
    const size_t prop_key = _data.porosity_key_index;
    std::vector<double> property(grid.n_active_cells());
    const std::string keyword = "PORO";
    size_t cell_index = 0;
    for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
      property[cell_index++] = _data.cell_properties[prop_key][cell->index()];

    IO::VTKWriter::add_data(property, keyword, out);
  }

  // save custom keywords
  for (std::size_t ivar=0; ivar<_data.output_flow_properties.size(); ++ivar)
  {
    const size_t prop_key = _data.output_flow_properties[ivar];
    std::vector<double> property(grid.n_active_cells());
    const std::string keyword = _data.property_names[prop_key];
    size_t cell_index = 0;
    for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
      property[cell_index++] = _data.cell_properties[prop_key][cell->index()];

    IO::VTKWriter::add_data(property, keyword, out);
  }

  // // save multiscale flow data
  // if (!data.ms_flow_data.partitioning.empty())
  // {
  //   IO::VTKWriter::add_data(data.ms_flow_data.partitioning, "partitioning-flow", out);
  //   saveMultiScaleSupport(data.ms_mech_data, grid.n_cells(), "support-flow-", out);
  // }
  // if (!data.ms_mech_data.partitioning.empty())
  // {
  //   const auto & ms = data.ms_mech_data;
  //   IO::VTKWriter::add_data(ms.partitioning, "partitioning-mech", out);
  //   IO::VTKWriter::enter_section_point_data(grid.n_vertices(), out);
  //   // saveMultiScaleSupport(msm, grid.n_vertices(), "support-mech-", out);

  //   // I'm just gonna save vertices and boundaries
  //   std::vector<int> output(grid.n_vertices(), 0);
  //   for (std::size_t coarse = 0; coarse < ms.n_coarse; coarse++)
  //   {
  //     for (const size_t i : ms.support_boundary[coarse])
  //       output[i] = 1;
  //   }
  //   for (const size_t i : ms.centroids)
  //     output[i] = 2;
  //   IO::VTKWriter::add_data(output, "support-mech", out);
  // }
  out.close();
}

void OutputDataVTK::save_reservoir_mechanics_data_(const std::string & fname) const
{
  logging::log() << "writing " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());
  // IO::VTKWriter::write_geometry(m_mech_grid, out);
  out << "# vtk DataFile Version 2.0 \n";
  out << "3D Grid\n";
  out << "ASCII \n \n";
  out << "DATASET UNSTRUCTURED_GRID \n";

  auto & grid = _data.geomechanics_grid;
  const size_t n_points = (_data.grid_vertices_after_face_split.empty()) ?
      grid.n_vertices() :
      _data.grid_vertices_after_face_split.size();

  out << "POINTS" << "\t" << n_points << " float" << std::endl;
  if (_data.grid_vertices_after_face_split.empty())
    for (const auto & vertex : grid.vertices())
      out << vertex << "\n";
  else
    for (const auto &vertex : _data.grid_vertices_after_face_split)
      out << vertex << "\n";

  const size_t n_entries_total = IO::VTKWriter::count_number_of_cell_entries_(grid);

  out << "CELLS " << grid.n_active_cells() << " " << n_entries_total << "\n";
  if (_data.grid_vertices_after_face_split.empty())
    for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
    {
      if ( cell->vtk_id() != angem::VTK_ID::GeneralPolyhedronID )
      {
        out << cell->n_vertices() << "\t";
        for (std::size_t i : cell->vertices())
          out << i << "\t";
        out << std::endl;
      }
      else
      {
        out << IO::VTKWriter::count_number_of_cell_entries_(*cell) << "\n";
        const auto & faces = cell->faces();
        out << faces.size() << "\n";
        for (const auto & face : faces)
        {
          const auto & vertices = face->vertices();
          out << vertices.size() << " ";
          for (const size_t v : vertices)
            out << v << " ";
          out << "\n";
        }
      }
    }
  else
    for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
    {
      auto cell_vertices = cell->vertices();
      auto true_vertices = _data.grid_cells_after_face_split[cell->index()];
      if ( cell->vtk_id() != angem::VTK_ID::GeneralPolyhedronID )
      {
        out << cell->n_vertices() << "\t";
        for (size_t i = 0; i < cell_vertices.size(); ++i)
          out << true_vertices[i] << "\t";
        out << "\n";
      }
      else
      {
        out << IO::VTKWriter::count_number_of_cell_entries_(*cell) << "\n";
        const auto & faces = cell->faces();
        out << faces.size() << "\n";
        for (const auto & face : faces)
        {
          const auto & vertices = face->vertices();
          out << vertices.size() << " ";
          for (const size_t v : vertices)
          {
            const size_t idx =  std::distance( cell_vertices.begin(),
                                  std::find( cell_vertices.begin(), cell_vertices.end(),
                                             v));
            out << true_vertices[idx] << " ";
          }
          out << "\n";
        }
      }
    }
  out << "CELL_TYPES" << "\t" << grid.n_active_cells() << "\n";
  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
    out << cell->vtk_id() << "\n";

  // save keywords
  // for (std::size_t ivar=0; ivar<_data.rockPropNames.size(); ++ivar)
  // {
  //   vector<double> property(grid.n_cells());
  //   const string keyword = data.rockPropNames[ivar];
  //   for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
  //     property[cell->index()] = data.vsCellRockProps[cell->index()].v_props[ivar];

  //   IO::VTKWriter::add_data(property, keyword, out);
  // }

  // // save multiscale flow data
  // if (!data.ms_flow_data.partitioning.empty())
  // {
  //   IO::VTKWriter::add_data(data.ms_flow_data.partitioning, "partitioning-flow", out);
  //   saveMultiScaleSupport(data.ms_mech_data, grid.n_cells(), "support-flow-", out);
  // }
  // if (!data.ms_mech_data.partitioning.empty())
  // {
  //   const auto & ms = data.ms_mech_data;
  //   IO::VTKWriter::add_data(ms.partitioning, "partitioning-mech", out);
  //   IO::VTKWriter::enter_section_point_data(grid.n_vertices(), out);
  //   // saveMultiScaleSupport(msm, grid.n_vertices(), "support-mech-", out);

  //   // I'm just gonna save vertices and boundaries
  //   std::vector<int> output(grid.n_vertices(), 0);
  //   for (std::size_t coarse = 0; coarse < ms.n_coarse; coarse++)
  //   {
  //     for (const size_t i : ms.support_boundary[coarse])
  //       output[i] = 1;
  //   }
  //   for (const size_t i : ms.centroids)
  //     output[i] = 2;
  //   IO::VTKWriter::add_data(output, "support-mech", out);
  // }
  out.close();
}



void OutputDataVTK::save_dfm_data(const std::string & fname) const
{
  logging::log() << "Saving DFM mesh file: " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());
  const auto & grid = _data.fracture_grid;
  IO::VTKWriter::write_surface_geometry(grid.get_vertices(), grid.get_polygons(), out);

  // save face markers
  IO::VTKWriter::enter_section_cell_data(grid.n_polygons(), out);
  {
    std::vector<double> property(grid.n_polygons());
    const std::string keyword = "Marker";
    for (auto face = grid.begin_polygons(); face != grid.end_polygons(); ++face)
      property[face.index()] = face.marker();
    IO::VTKWriter::add_data(property, keyword, out);
  }

  out.close();
}


void OutputDataVTK::save_edfm_data(const std::string & fname) const
{
  // std::cout << "Saving EDFM mesh file: " << fname << std::endl;
  // std::ofstream out;
  // out.open(fname.c_str());

  // // write vtk data
  // std::size_t n_efrac_vertices = 0;
  // // make up a vector of all sda vertices
  // for (const auto & efrac : data.vEfrac)
  //   n_efrac_vertices += efrac.mesh.n_vertices();

  // std::vector<angem::Point<3,double>> efrac_verts(n_efrac_vertices);
  // std::size_t ivertex = 0;
  // for (const auto & efrac : data.vEfrac)
  //   for (const auto & p : efrac.mesh.get_vertices())
  //   {
  //     efrac_verts[ivertex] = p;
  //     ivertex++;
  //   }

  // std::size_t n_efrac_elements = 0;
  // std::size_t vind_size_total = 0;

  // for (const auto & efrac : data.vEfrac)
  // {
  //   n_efrac_elements += efrac.mesh.n_polygons();
  //   for (const auto & poly : efrac.mesh.get_polygons())
  //     vind_size_total += poly.size();
  // }

  // std::vector<std::vector<std::size_t>> efrac_cells(n_efrac_elements);
  // std::size_t ielement = 0;
  // std::size_t shift = 0;
  // for (const auto & efrac : data.vEfrac)
  // {
  //   for (const auto & cell : efrac.mesh.get_polygons())
  //   {
  //     efrac_cells[ielement].resize(cell.size());
  //     for (short v=0; v<cell.size(); ++v)
  //       efrac_cells[ielement][v] = shift + cell[v];
  //     ielement++;
  //   }
  //   shift += efrac.mesh.n_vertices();
  // }

  // IO::VTKWriter::write_surface_geometry(efrac_verts, efrac_cells, out);
  // IO::VTKWriter::write_surface_geometry(_data.edfm_grid.get_vertices(),
  //                                       _data.edfm_grid.get_polygons(), out);

  // // save face markers
  // const auto & grid = _data.edfm_grid;
  // IO::VTKWriter::enter_section_cell_data(grid.n_polygons(), out);
  // {
  //   std::vector<double> property(grid.n_polygons());
  //   const std::string keyword = "Marker";
  //   for (auto face = grid.begin_polygons(); face != grid.end_polygons(); ++face)
  //     property[face.index()] = face.marker();
  //   IO::VTKWriter::add_data(property, keyword, out);
  // }

  // out.close();
}


// void OutputDataVTK::saveMultiScaleSupport(const multiscale::MultiScaleOutputData & ms,
//                                           const std::size_t                        size,
//                                           const std::string                      & prefix,
//                                           std::ofstream                          & out)
// {
//   // coarse nodes
//   std::vector<int> support_value(size);
//   for (std::size_t coarse = 0; coarse < ms.n_coarse; coarse++)
//   {
//     for (std::size_t i=0; i<size; ++i)
//     {
//       if (i == ms.centroids[coarse])
//         support_value[i] = 3;
//       else if (ms.support_boundary[coarse].find(i) != ms.support_boundary[coarse].end())
//         support_value[i] = 2;
//       else if (ms.support_internal[coarse].find(i) != ms.support_internal[coarse].end())
//         support_value[i] = 1;
//       else
//         support_value[i] = 0;
//     }

//     IO::VTKWriter::add_data(support_value, prefix + std::to_string(coarse), out);
//   }  // end coarse loop
// }


void OutputDataVTK::save_wells_(const std::string & fname) const
{
  logging::log() << "writing " << fname << std::endl;
  IO::VTKWriter::write_well_trajectory(_data.well_vertices.points, _data.well_vertex_indices, fname);
}

}
