#include "OutputDataGPRS.hpp"
#include "VTKWriter.hpp"     // provides IO::VTKWriter
#include "logger/Logger.hpp"  // provides logging::log


namespace gprs_data
{

const int n_entries_per_line = 10;
static constexpr double transmissibility_conversion_factor = 0.0085267146719160104986876640419948;


OutputDataGPRS::OutputDataGPRS(const SimData & data, const GPRSOutputConfig config)
    :
    _data(data),
    _config(config)
{}


void OutputDataGPRS::write_output(const std::string & output_path) const
{
  _output_path = output_path;
  // save flow data
  save_flow_data_(_output_path + "/" + _config.flow_cv_file,
                  _output_path + "/" + _config.flow_connection_file);

  if (_data.has_mechanics)
    save_geomechanics_data_();

  // std::cout << "save custom keyswords" << std::endl;
  // saveGeomechDataNewKeywords(output_path + data.config.domain_file);

  if (!_data.sda_data.empty())
    save_embedded_fractures_(_output_path + "/" + _config.efrac_file);

  // std::cout << "save mech boundary conditions: "
  //           << output_path + data.config.bcond_file
  //           << std::endl;
  // saveBoundaryConditions(output_path + data.config.bcond_file);

  if (_data.dfm_faces.size() > 0 && _data.has_mechanics)
    save_discrete_fracture_properties_(output_path + "/" + _config.discrete_frac_file);

  if (!_data.wells.empty())
    saveWells(output_path + "/" + _config.wells_file);

  // // multiscale
  // if (data.ms_flow_data.partitioning.size() > 0)
  //   saveFlowMultiScaleData(output_path + data.config.flow_ms_file);
  // if (data.ms_mech_data.partitioning.size() > 0)
  //   saveMechMultiScaleData(output_path + data.config.mech_ms_file);
}

void OutputDataGPRS::save_flow_data_(const std::string cv_file, const std::string con_file) const
{
  {  // Write cell data
    std::ofstream out;
    out.open(cv_file.c_str());
    logging::log() << "saving flow Control Volume data: " << cv_file << std::endl;
    save_control_volume_data_(out);
    out.close();
  }
  {  // write face data
    std::ofstream out;
    out.open(con_file.c_str());
    logging::log() << "saving flow Connection data: " << con_file << std::endl;
    save_trans_data_(out);
    save_trans_update_formulas_(out);
    // write transmissibility update formulas
    out.close();
  }
}

void OutputDataGPRS::save_control_volume_data_(std::ofstream & out) const
{
  const auto & cvs = _data.cv_data;
  ///// OUTPUT Dimensions /////
  out << "DIMENS" << std::endl;
  out << cvs.size() << "\t" << 1 << "\t" << 1 << "\t" << std::endl;
  out << "/" << std::endl << std::endl;

  ///// OUTPUT Volumes /////
  out << "VOLUME" << std::endl;
  for (const auto & cv : cvs)
    out << cv.volume << std::endl;
  out << "/" << std::endl << std::endl;

  ///// OUTPUT Porosity /////
  out << "PORO" << std::endl;
  for (const auto & cv : cvs)
    out << cv.porosity << std::endl;
  out << "/" << std::endl << std::endl;

  ///// OUTPUT Depth  /////
  out << "DEPTH" << std::endl;
  for (const auto & cv : cvs)
    out << -cv.center(2) << std::endl;
  out << "/" << std::endl << std::endl;

  // additional data (if any)
  for (std::size_t ivar=0; ivar<_data.output_flow_properties.size(); ++ivar)
  {
    const size_t prop_key = _data.output_flow_properties[ivar];
    const std::string keyword = _data.property_names[prop_key];
    out << keyword << std::endl;
    for (const auto & cv : cvs)
        out << cv.custom[ivar] << std::endl;
    out << "/" << std::endl << std::endl;
  }
}

void OutputDataGPRS::save_trans_data_(std::ofstream & out) const
{
  const auto & cons = _data.flow_connection_data;

  /* OUTPUT Transmissibility */
  out << "TPFACONNS" << std::endl;
  // n nonzero connection
  std::size_t n_connections = 0;
  for (const auto & con : cons)
  {
    const double transissibility = std::fabs(con.coefficients[0]) * transmissibility_conversion_factor;
    if (transissibility > 1e-10) n_connections++;
  }
  out << n_connections << std::endl;
  for (const auto & con : cons)
  {
    assert( con.elements.size() == 2 );
    assert( con.coefficients.size() == 2 );

    const double transissibility = std::fabs(con.coefficients[0]) * transmissibility_conversion_factor;
    if (transissibility > 1e-10)
      out << con.elements[0] << "\t" << con.elements[1] << "\t"
          << std::scientific << transissibility << std::defaultfloat << std::endl;
  }
  out << "/" << std::endl;
}

void OutputDataGPRS::save_trans_update_formulas_(std::ofstream & out) const
{
  out << "GMUPDATETRANS\n";
  const auto & cons = _data.flow_connection_data;

  for (size_t icon = 0; icon < cons.size(); ++icon)
  {
    const auto & con = cons[icon];
    out << icon << "\t";
    switch (con.type)
    {
      case discretization::ConnectionType::matrix_matrix:
        out << con.type << "\t"
            << con.elements[0] << "\t" << con.update_formula[0] << "\t"
            << con.elements[1] << "\t" << con.update_formula[1] << "\n";
        break;
      case discretization::ConnectionType::matrix_fracture:
        if (con.update_formula.size() < 4)
          throw std::runtime_error("FIXME: invalid formula for trans update in M-F case");
        out << con.type << "\t"
            << con.elements[0] << "\t" << con.update_formula[0]
            << "\t" << con.update_formula[1] << "\t" << con.update_formula[2]
            << "\t" << con.elements[1] << "\t" << con.update_formula[3] << "\n";
        break;
      case discretization::ConnectionType::fracture_fracture:
        out << con.type <<"\t"
            << con.elements[0] << "\t"
            << con.update_formula[0] <<"\t"
            << con.update_formula[1] << "\t"
            << con.update_formula[2] << "\t"
            << con.elements[1] << "\t" << con.update_formula[3]
            << "\t" << con.update_formula[4]
            << "\t" << con.update_formula[5] << "\t";
        // remaining elements in the connection
        out << con.all_elements.size() << "\t";
        size_t shift = 6;
        for (size_t j = 0; j < con.all_elements.size(); ++j)
        {
          out << con.all_elements[j] << "\t";
          out << con.update_formula[shift++] << " ";  // T_j / K_j
          out << con.update_formula[shift++] << " ";  // volume factor
          out << con.update_formula[shift++] << "\t";  // K_j
        }
        out << "\n";
        break;
    }
  }
  out << "/\n";
}

void OutputDataGPRS::save_geometry_() const
{
  const std::string outstring = _output_path + "/" + _config.geometry_file;
  logging::log() << "saving mechaanics geometry: " << outstring << std::endl;

  // gprs output
  std::ofstream out;
  out.open(outstring.c_str());
  out << "GMDIMS" << "\n";

  const auto & grid = _data.geomechanics_grid;

  // number of vertices: either take from grid or from split constainer
  size_t nv = grid.n_vertices();
  if (!_data.grid_vertices_after_face_split.empty())
    nv = _data.grid_vertices_after_face_split.size();

  size_t n_active_faces = 0;
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
    n_active_faces++;

  out << nv << "\t"
      << grid.n_active_cells() << "\t"
      << n_active_faces;
  out << "/" << "\n\n";

  // write vertex coordinates
  out.precision(6);
  out << "GMNODE_COORDS" << "\n";
  if (_data.grid_vertices_after_face_split.empty())
  {
    for (const auto & vertex : grid.vertices())
      out << vertex[0] << "\t" << vertex[1] << "\t" << vertex[2] << "\n";
  }
  else
  {
    for (const auto &vertex : _data.grid_vertices_after_face_split)
      out << vertex[0] << "\t" << vertex[1] << "\t" << vertex[2] << "\n";
  }
  out << "/" << std::endl << std::endl;

  // Elemenent data
  save_cell_geometry_(out, grid);

  // geomechanics -> flow cell mapping
  out << "GMCELL_TO_FLOWCELLS" << "\n";
  for (const auto & flow_cells : _data.gmcell_to_flowcells)
  {
    const std::size_t n_connected_elements = flow_cells.size();
    if (n_connected_elements == 0)
      out << 1 << "\t" << -1 << std::endl;
    else
    {
      out << n_connected_elements << "\t";
      for (const std::size_t ielement : flow_cells)
        out << ielement + 1 << "\t";
      out << std::endl;
    }
  }
  out << "/\n\n";

  // Face data
  save_face_geometry_(out, grid);

  // map faces to cells
  out << "GMFACE_GMCELLS" << std::endl;
  for (auto face=grid.begin_active_faces(); face!=grid.end_active_faces(); ++face)
  {
    const std::vector<const mesh::Cell*> neighbors = face->neighbors();
    out << neighbors.size() << "\t";
    for (const mesh::Cell* neighbor : neighbors)
      out << _data.mech_numbering->cell_dof(neighbor->index()) + 1 << "\t";
    out << "\n";
  }
  out << "/\n\n";
}


void OutputDataGPRS::save_geomechanics_keywords_() const
{
  // write domain properties
  std::ofstream out;
  const std::string file_name = _output_path + "/" + _config.mechanics_kwd_file;
  logging::log() << "saving geomechanics keywords: " << file_name << std::endl;
  out.open(file_name.c_str());

  const auto & grid = _data.geomechanics_grid;
  for (std::size_t ivar=0; ivar<_data.output_mech_properties.size(); ++ivar)
  {
    const size_t prop_key = _data.output_mech_properties[ivar];
    const std::string keyword = _data.property_names[prop_key];
    out << keyword << "\n";
    size_t cnt = 0;
    for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
    {
      const std::size_t icell = cell->index();
      out << _data.cell_properties[prop_key][cell->index()] << "\t";
      if (++cnt % n_entries_per_line == 0)
        out << "\n";
    }
    out << "/\n\n";
  }

  out.close();
}


void OutputDataGPRS::save_embedded_fractures_(const std::string file_name) const
{
  logging::log() << "saving mech embedded fracture data: " << file_name << std::endl;
  std::ofstream out;
  out.open(file_name.c_str());

  // write SDA Geomechanics cells
  out << "GM_EFRAC_CELLS" << std::endl;
  for (const auto & frac : _data.sda_data)
  {
    out << frac.cells.size() << std::endl << "\t";
    for (size_t i=0; i < frac.cells.size(); ++i)
    {
      const double mech_cell = _data.mech_numbering->cell_dof(frac.cells[i]);
      out << mech_cell + 1 << "\t";
      if ((i+1) % n_entries_per_line == 0)
        out << "\n";
      if (i == frac.cells.size() - 1)
        out << "\n";
    }
  }
  out << "/\n\n";

  // write SDA Flow cells for each geomechanics cell in embedded fractures defined in
  // "GM_EFRAC_CELLS".
  // The order is the same as in the order of geomechanics cells in "GM_EFRAC_CELLS".
  out << "GM_EFRAC_TO_FLOWCELLS" << std::endl;
  for (const auto & frac : _data.sda_data)
  {
    for (size_t i=0; i < frac.cells.size(); ++i)
    {
      const size_t mech_cell = _data.mech_numbering->cell_dof(frac.cells[i]);
      std::vector<size_t> flow_cells = _data.gmcell_to_SDA_flowcells[mech_cell]; // flow cells of geo cells in GM_EFRAC_CELLS
      if (flow_cells.size() == 0) // uncoupled case
        out << 0 << "\t" << -1 << std::endl;
      else // coupled case
      {
        out << flow_cells.size() << "\t";
        for (const std::size_t ielement : flow_cells)
          out << ielement + 1 << "\t";
        out << std::endl;
      }
    }
  }
  out << "/\n\n";

  // Define a fracture region for each fracture.
  // The fractures in the same fracture region is allocated with the
  // same tabulated data for fracture permeability and volume multiplier
  // defined with "GMCONTACT_NORMAL_PROPS"
  out << "GM_EFRAC_REGION" << std::endl;
  for (const auto & frac : _data.sda_data)
    out << frac.region + 1 << "\n";
  out << "/\n\n";

  // coordinates of a point in frac plane for each SDA cell
  out << "GM_EFRAC_POINTS" << "\n";
  for (const auto & frac : _data.sda_data)
      for (const auto & point : frac.points)
      out << point << "\n";
  out << "/\n\n";

  // sda fracture dip in each fracture cell
  out << "GM_EFRAC_DIP" << "\n";
  for (const auto & frac : _data.sda_data)
    for (size_t i=0; i < frac.dip.size(); ++i)
    {
      out << frac.dip[i] << "\t";
      if ((i+1) % n_entries_per_line == 0)
        out << "\n";
    }
  out << "/\n\n";

  out << "GM_EFRAC_STRIKE" << "\n";
  for (const auto & frac : _data.sda_data)
    for (size_t i=0; i < frac.strike.size(); ++i)
    {
      out << frac.strike[i] << "\t";
      if ((i+1)%n_entries_per_line == 0) out << "\n";
    }
  out << "/\n\n";

  out << "GM_EFRAC_COHESION" << "\n";
  for (const auto & frac : _data.sda_data)
    out << frac.cohesion << "\n";
  out << "/\n\n";

  out << "GM_EFRAC_FRICTION" << "\n";
  for (const auto & frac : _data.sda_data)
    out << frac.friction_angle << "\n";
  out << "/\n\n";

  out << "GM_EFRAC_DILATION" << "\n";
  for (const auto & frac : _data.sda_data)
    out << frac.dilation_angle << "\n";
  out << "/\n\n";

  // fracture conductivity of each fracture (md-m)
  out << "GM_EFRAC_CONDUCTIVITY" << std::endl;
  for (const auto & frac : _data.sda_data)
    out << frac.conductivity << "\n";
  out << "/\n\n";

  out.close();
}


void OutputDataGPRS::save_geomechanics_boundary_conditions_() const
{
  const std::string file_name = _output_path + "/" + _config.bcond_file;
  logging::log() << "saving geomech boundary conditions: " << file_name << std::endl;

  std::ofstream out(file_name);
  // save dirichlet BC's
  if ( !_data.dirichlet_indices[0].empty() )
    save_dirichlet_component_vertices(0, "X", out);
  if ( !_data.dirichlet_indices[1].empty() )
    save_dirichlet_component_vertices(1, "Y", out);
  if ( !_data.dirichlet_indices[2].empty() )
    save_dirichlet_component_vertices(2, "Z", out);

  // save dirichlet constraints
  if (!_data.boundary_constraints.empty())
  {
    out << "GMEQCONSTRAINT" << "\n";
    for (size_t igroup = 0; igroup < _data.boundary_constraints.size(); ++igroup)
    {
      for (size_t i = 1; i < _data.boundary_constraints[igroup].nodes.size(); ++i)
      {
        out << "2 ";
        out << _data.boundary_constraints[igroup].nodes[i-1] + 1 << " "
            << _data.boundary_constraints[igroup].nodes[i]   + 1 << " ";
        out << _data.boundary_constraints[igroup].components[i-1] + 1 << " "
            << _data.boundary_constraints[igroup].components[i] + 1  << " ";
        out << 1 << " " << -1 << " ";
      out << _data.boundary_constraints[igroup].penalty << "\n";
      }
    }
    out << "/\n\n";
  }

  if ( !_data.neumann_face_indices.empty() )
  {
    out << "GMFACE_TRACTION_TXYZ" << "\n";

    for (std::size_t i=0; i<_data.neumann_face_indices.size(); ++i)
    {
      out << _data.mech_numbering->face_dof(_data.neumann_face_indices[i]) + 1 << "\t";
      out << _data.neumann_face_traction[i] << "\n";
    }
    out << "/\n\n";
  }

  out.close();
}

void OutputDataGPRS::save_dirichlet_component_vertices(const size_t comp,
                                                       const std::string comp_name,
                                                       std::ofstream & out) const
{
  out << "GMNODE_BCDISP" << comp_name << "\n";
  for (std::size_t i=0; i<_data.dirichlet_indices[comp].size(); ++i)
  {
    out << _data.dirichlet_indices[comp][i] + 1 << "\t";
    out << _data.dirichlet_values[comp][i] << "\n";
  }
  out << "/\n\n";
}

void OutputDataGPRS::save_discrete_fracture_properties_(const std::string file_name) const
{
  std::ofstream out;
  out.open(file_name.c_str());
  logging::log() << "saving geomech DFM file: " << file_name << std::endl;

  logging::debug() << "write all fractured faces\n";
  out << "GMFACE_FRACTURE_TO_FLOWCELL" << std::endl;
  for (const auto & it : _data.dfm_faces)
  {
    out << _data.mech_numbering->face_dof(it.first) + 1 << "\t";
    if (it.second.coupled)
      out << _data.flow_numbering->face_dof(it.first) + 1 << std::endl;
    else
      out << -1 << std::endl;
  }
  out << "/\n\n";

  out << "GMFACE_FRACTURE_CONDUCTIVITY" << std::endl;
  for (const auto facet_it : _data.dfm_faces)
    out << facet_it.second.conductivity << std::endl;
  out << "/\n\n";

  out << "GMFACE_FRACTURE_REGION" << std::endl;
  for (const auto & it : _data.dfm_faces)
    out << it.second.region + 1 << "\n";
  out << "/\n\n";

  out << "GMFACE_FRACTURE_GROUP" << std::endl;
  for (const auto & it : _data.dfm_faces)
    out << 1 << "\n";
  out << "/\n\n";

  out.close();
}


void OutputDataGPRS::saveWells(const std::string file_name) const
{
  logging::log() << "Save wells: " << file_name << std::endl;
  std::ofstream out;
  out.open(file_name.c_str());

  out << "WELSPECS" << std::endl;
  for (const auto & well : _data.wells)
  {
    out << well.name << "\tGROUP1\t";
    assert( !well.segment_data.empty() && "Well does not contain connection data" );
    out << well.segment_data.front().dof << " 1 ";
    // reference depth
    out << -well.reference_depth << " /" << std::endl;
  }
  out << "/" << std::endl << std::endl;

  out << "COMPDAT" << std::endl;
  out << "-- name\tcell\tidk\tidk\tidk\topen\tidk\tWI\trad" << std::endl;
  for (const auto & well : _data.wells)
  {
    // for (std::size_t i=0; i<well.connected_volumes.size(); ++i)
    for (const auto & segment : well.segment_data)
    {
      out << well.name << "\t";
      out << segment.dof + 1 << "\t";
      // j, k1:k2 open sat_table_number
      out << "1\t1\t1\tOPEN\t1*\t";
      out << segment.wi * transmissibility_conversion_factor << "\t";
      out << 2*well.radius << "\t";  // adgprs needs well diameter
      out << "/" << std::endl;
    }
  }
  out << "/" << std::endl << std::endl;

  out.close();
}


void OutputDataGPRS::saveFlowMultiScaleData(const std::string file_name)
{

}


void OutputDataGPRS::saveMechMultiScaleData(const std::string file_name)
{
  // std::ofstream out;
  // out.open(file_name.c_str());
  // const auto & ms = data.ms_mech_data;

  // // save partitioing
  // out << "GMMSPARTITIONING";
  // for (std::size_t i=0; i < ms.partitioning.size(); ++i)
  // {
  //   if (i % n_entries_per_line == 0) out << endl;
  //   out << ms.partitioning[i] << " ";
  // }
  // out << "/" << endl << endl;

  // // save support
  // out << "GMMSSUPPORT ";
  // for (std::size_t i=0; i < ms.n_coarse; ++i)
  // {
  //   out << endl;
  //   out << ms.support_internal[i].size() << " "  // number of cells (centroid)
  //       << ms.support_boundary[i].size() << " "; // number of boundary nodes

  //   // first write the centroid (vertex)
  //   out << ms.centroids[i] << " ";

  //   // internal cells
  //   size_t counter = 3;
  //   for (const size_t cell : ms.support_internal[i])
  //   {
  //     if (counter++ % n_entries_per_line == 0) out << endl;
  //     out << cell << " ";
  //   }

  //   // boundary nodes
  //   for (const size_t vertex : ms.support_boundary[i])
  //   {
  //     if (counter++ % n_entries_per_line == 0) out << endl;
  //     out << vertex << " ";
  //   }
  // }

  // out << "/" << endl << endl;

  // out.close();
}

void OutputDataGPRS::save_geomechanics_data_() const
{
  save_geometry_();
  // save computed element data: shape functions and gradients, gauss weights, JxW
  save_fem_data_();

  save_geomechanics_keywords_();

  save_geomechanics_boundary_conditions_();
}

void save_cell_vertices(std::ofstream & out, const std::vector<size_t> &vertices, const int vtk_id)
{
  out << vertices.size() << "\t";
  switch (vtk_id)
  {
    case 25: // super wierd element 25
      {
        for (int j = 0; j < 8; j++)
          out << vertices[j] + 1 << "\t";
        out << vertices[8] + 1 << "\t";
        out << vertices[11] + 1 << "\t";
        out << vertices[13] + 1 << "\t";
        out << vertices[9] + 1 << "\t";

        out << vertices[16] + 1 << "\t";
        out << vertices[18] + 1 << "\t";
        out << vertices[19] + 1 << "\t";
        out << vertices[17] + 1 << "\t";

        out << vertices[10] + 1 << "\t";
        out << vertices[12] + 1 << "\t";
        out << vertices[14] + 1 << "\t";
        out << vertices[15] + 1 << "\t";
        break;
      }
    case 26:
      {
        for (int j = 0; j < 6; j++)
          out << vertices[j] + 1 << "\t";

        out << vertices[6] + 1 << "\t";
        out << vertices[9] + 1 << "\t";
        out << vertices[7] + 1 << "\t";

        out << vertices[12] + 1 << "\t";
        out << vertices[14] + 1 << "\t";
        out << vertices[13] + 1 << "\t";

        out << vertices[8] + 1 << "\t";
        out << vertices[10] + 1 << "\t";
        out << vertices[11] + 1 << "\t";
        break;
      }
    default:
      {
        for (const auto vertex : vertices)
          out << vertex + 1 << "\t";
        break;
      }
  }
  out << "\n";
}

void OutputDataGPRS::save_cell_geometry_(std::ofstream & out, const mesh::Mesh & grid) const
{
  out << "GMCELL_NODES" << "\n";
  if (_data.grid_cells_after_face_split.empty())
    for (auto cell=grid.begin_active_cells(); cell!=grid.end_active_cells(); ++cell)
      save_cell_vertices(out, cell->vertices(), cell->vtk_id());
  else
    for (auto cell=grid.begin_active_cells(); cell!=grid.end_active_cells(); ++cell)
      save_cell_vertices(out, _data.grid_cells_after_face_split[cell->index()], cell->vtk_id());

  out << "/" << "\n\n";

  // Cell types //
  out << "GMCELL_TYPE" << "\n";
  for (auto cell=grid.begin_active_cells(); cell!=grid.end_active_cells(); ++cell)
    out <<  cell->vtk_id() << "\n";
  out << "/" << "\n\n";
}

void export_point_grads(std::ofstream & out, const discretization::FiniteElementData & data)
{
  for (const auto & point_data : data.points)  // loop q points
  {
    for (const auto & grad : point_data.grads)  // loop vertices
      out << grad << "\t";                      // writes 3 values (x,y,z)
    out << "\n";
  }
 
}

void OutputDataGPRS::save_face_geometry_(std::ofstream & out, const mesh::Mesh & grid) const
{
  out << "GMFACE_NODES\n";
  for (auto face=grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
  {
    const std::vector<std::size_t> ivertices = face->vertices();
    out << ivertices.size() << "\t";
    for (const std::size_t ivertex : ivertices)
      out << ivertex + 1 << "\t";
    out << "\n";
  }
  out << "/" << "\n\n";

  out << "GMFACE_TYPE" << std::endl;
  for (auto face = grid.begin_active_faces(); face!=grid.end_active_faces(); ++face)
    out << face->vtk_id() << std::endl;
  out << "/" << "\n\n";
}

void OutputDataGPRS::save_fem_data_() const
{
  if (_data.fe_cell_data.empty() && _data.fe_face_data.empty())
    return;

  const std::string file_name = _output_path + "/" + _config.fem_file;
  logging::log() << "saving FEM data: " << file_name << std::endl;
  std::ofstream out;
  out.open(file_name.c_str());
  // save cell data
  const std::vector<discretization::FiniteElementData> & cells = _data.fe_cell_data;
  out << "GMCELL_GAUSS_WEIGHTS" << "\n";
  for (const auto & cell : cells)
    if (!cell.points.empty())
  {
    out << _data.mech_numbering->cell_dof(cell.element_index) + 1 << "\t";
    out << cell.points.size() << "\t";
    for (const auto & point : cell.points)
      out << point.weight << "\t";
    out << "\n";
  }
  out << "/\n\n";
  // cell center
  out << "GMCELL_GAUSS_WEIGHTS_CENTER" << "\n";
  for (const auto & cell : cells)
    if (!cell.points.empty())
    {
      out << cell.center.weight << "\t";
      out << "\n";
    }
  out << "/\n\n";

  out << "GMCELL_SHAPE_VALUES" << "\n";
  for (const auto & cell : cells)
    if (!cell.points.empty())
    {
      for (const auto &point : cell.points)
      {
        for (const double value : point.values)
          out << value << "\t";
        out << "\n";
      }
    }
  out << "/\n\n";
  // center
  out << "GMCELL_SHAPE_VALUES_CENTER" << "\n";
  for (const auto & cell : cells)
    if (!cell.points.empty())
    {
      for (const double value : cell.center.values)
          out << value << "\t";
      out << "\n";
    }
  out << "/\n\n";

  out << "GMCELL_SHAPE_GRADS" << "\n";
  for (const auto & cell : cells)
    if (!cell.points.empty())
    {
      export_point_grads(out, cell);
    }
  out << "/\n\n";
  // center
  out << "GMCELL_SHAPE_GRADS_CENTER" << "\n";
  for (const auto & cell : cells)
    if (!cell.points.empty())
    {
      for (const angem::Point<3,double> & grad : cell.center.grads)
          out << grad << "\t";
      out << "\n";
    }
  out << "/\n\n";

  // SAVE FACE DATA
  const auto & faces = _data.fe_face_data;  // fe values and gradients for grid faces
  if (faces.empty())
  {
    out.close();
    return;
  }

  // face gauss weights
  out << "GMFACE_GAUSS_WEIGHTS" << "\n";
  for (const auto & face : faces)
    if (!face.points.empty())
    {
      out << _data.mech_numbering->face_dof(face.element_index) + 1 << "\t";
      out << face.points.size() << "\t";
      for (const auto & point : face.points)
        out << point.weight << "\t";
      out << "\n";
    }
  out << "/\n\n";

  // center detJ
  out << "GMFACE_GAUSS_WEIGHTS_CENTER" << "\n";
  for (const auto & face : faces)
    if (!face.points.empty())
      out << face.center.weight << "\n";
  out << "/\n\n";

  // face shape function values
  out << "GMFACE_SHAPE_VALUES" << "\n";
  for (const auto & face : faces)
    if (!face.points.empty())
    {
      for (const auto &point : face.points)
      {
        for (const double value : point.values)
          out << value << "\t";
        out << "\n";
      }
    }
  out << "/\n\n";

  out << "GMFACE_SHAPE_VALUES_CENTER" << "\n";
  for (const auto & face : faces)
    if (!face.points.empty())
    {
      for (const double value : face.center.values)
        out << value << "\t";
      out << "\n";
    }
  out << "/\n\n";

  // output face shape grads
  out << "GMFACE_SHAPE_GRADS" << "\n";
  for (const auto & face : faces)
    if (!face.points.empty())
      export_point_grads(out, face);
  out << "/\n\n";

  out << "GMFACE_SHAPE_GRADS_CENTER" << "\n";
  for (const auto & face : faces)
    if (!face.points.empty())
    {
      for (const angem::Point<3,double> & grad : face.center.grads)
        out << grad << "\t";
      out << "\n";
    }
  out << "/\n\n";

  out << "GMFRACTURE_SHAPE_GRADS" << "\n";
  for (size_t iface = 0; iface < _data.fe_frac_data.size(); ++iface)
    if ( !_data.fe_frac_data[iface].empty() )
    {
      const auto face = _data.geomechanics_grid.face(iface);
      const auto neighbors = face.neighbors();
      auto & data = _data.fe_frac_data[iface];
      if (data[0].element_index == neighbors[0]->index())
      {
        assert( data[0].points.size() == data[1].points.size() );
        export_point_grads(out, data[0]);
        export_point_grads(out, data[1]);
      }
      else
      {
        if ( data[1].element_index != neighbors[0]->index() )
          throw std::runtime_error("something went wront");
        export_point_grads(out, data[1]);
        export_point_grads(out, data[0]);
      }
    }

  out << "/\n\n";

  out.close();
}

}
