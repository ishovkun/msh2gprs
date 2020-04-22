#include "OutputDataGPRS.hpp"
#include "VTKWriter.hpp"

#include <sys/stat.h>

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

  save_geomechanics_data_();

  // std::cout << "save custom keyswords" << std::endl;
  // saveGeomechDataNewKeywords(output_path + data.config.domain_file);

  // if (!data.vEfrac.empty())
  // {
  //   std::cout << "save embedded fractures" << std::endl;
  //   saveEmbeddedFractureProperties(output_path + data.config.efrac_file);
  // }

  // std::cout << "save mech boundary conditions: "
  //           << output_path + data.config.bcond_file
  //           << std::endl;
  // saveBoundaryConditions(output_path + data.config.bcond_file);

  // if (data.dfm_faces.size() > 0)
  // {
  //   std::cout << "save discrete fractures" << std::endl;
  //   saveDiscreteFractureProperties(output_path + data.config.discrete_frac_file);
  // }

  if (!_data.wells.empty())
  {
    std::cout << "save wells" << std::endl;
    saveWells(output_path + "/" + _config.wells_file);
  }

  // // flow discretization
  // std::cout << "save flow discretization" << std::endl;
  // flow::CalcTranses::save_output(data.flow_data, output_path);

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
    std::cout << "saving " << cv_file << std::endl;

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

    out.close();
  }
  {  // write face data
    std::ofstream out;
    out.open(con_file.c_str());
    std::cout << "saving " << con_file << std::endl;

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

    out.close();
  }
}

void OutputDataGPRS::save_geometry_() const
{
  const std::string outstring = _output_path + "/" + _config.geometry_file;
  std::cout << "saving " << outstring << std::endl;

  // gprs output
  std::ofstream out;
  out.open(outstring.c_str());
  out << "GMDIMS" << "\n";

  const auto & grid = _data.geomechanics_grid;
  size_t n_active_faces = 0;
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
    n_active_faces++;

  out << grid.n_vertices() << "\t"
              << grid.n_active_cells() << "\t"
              << n_active_faces;
  out << "/" << "\n\n";

  // write vertex coordinates
  out.precision(6);
  // std::cout << "write all coordinates\n";
  out << "GMNODE_COORDS" << "\n";
  for (const auto & vertex : grid.vertices())
      out << vertex[0] << "\t"
                  << vertex[1] << "\t"
                  << vertex[2] << "\n";
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
  //   if (data.is_fracture(face.marker()))  // timur want to retain neighbors of master frac face
  //   {
  //     if (face.index() == face.master_index())
  //     {
  //       const auto it_frac_face = data.dfm_faces.find(face.master_index());
  //       if (it_frac_face == data.dfm_faces.end())
  //       {
  //         throw std::runtime_error("bug in dfm connections");
  //       }
  //       const auto & neighbors = it_frac_face->second.neighbor_cells;

  //       out << neighbors.size() << "\t";
  //       for (const auto & neighbor : neighbors)
  //         out << neighbor + 1 << "\t";
  //       out << std::endl;
  //     }
  //   }
  //   else
  //   {
      const std::vector<const mesh::Cell*> neighbors = face->neighbors();
      out << neighbors.size() << "\t";
      for (const mesh::Cell* neighbor : neighbors)
        out << _data.mech_cell_numbering->cell_dof(neighbor->index()) + 1 << "\t";
      out << "\n";
  //   }
  }
  out << "/\n\n";

}


void OutputDataGPRS::save_geomechanics_keywords_() const
{
  // write domain properties
  std::ofstream out;
  const std::string file_name = _output_path + "/" + _config.mechanics_kwd_file;
  std::cout << "saving " << file_name << std::endl;
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


void OutputDataGPRS::saveEmbeddedFractureProperties(const std::string file_name)
{
  // std::cout  << "Writing SDA props" << std::endl;
  // std::ofstream geomechfile;
  // geomechfile.open(file_name.c_str());

  // geomechfile << "GM_EFRAC_CELLS" << std::endl;
  // for (const auto & efrac : data.vEfrac)
  // {
  //   geomechfile << efrac.cells.size() << std::endl << "\t";
  //   for (std::size_t i=0; i<efrac.cells.size(); ++i)
  //   {
  //     geomechfile << efrac.cells[i] + 1 << "\t";
  //     if ((i+1) % n_entries_per_line == 0)
  //       geomechfile << std::endl;
  //     if (i == efrac.cells.size() - 1)
  //       geomechfile << std::endl;
  //   }
  // }
  // geomechfile << "/" << std::endl << std::endl;

  // geomechfile << "GM_EFRAC_POINTS" << std::endl;
  // for (const auto & efrac : data.vEfrac)
  //   for (std::size_t i=0; i<efrac.points.size(); ++i)
  //     geomechfile << efrac.points[i] << std::endl;
  // geomechfile << "/" << std::endl << std::endl;

  // geomechfile << "GM_EFRAC_DIP" << std::endl;
  // for (const auto & efrac : data.vEfrac)
  //   for (std::size_t i=0; i<efrac.points.size(); ++i)
  //   {
  //     geomechfile << efrac.dip[i] << "\t";
  //     if ((i+1) % n_entries_per_line == 0)
  //       geomechfile << std::endl;
  //   }
  // geomechfile << "/" << std::endl << std::endl;

  // geomechfile << "GM_EFRAC_STRIKE" << std::endl;
  // for (const auto & efrac : data.vEfrac)
  //   for (std::size_t i=0; i<efrac.points.size(); ++i)
  //   {
  //     geomechfile << efrac.strike[i] << "\t";
  //     if ((i+1)%n_entries_per_line == 0) geomechfile << std::endl;
  //   }
  // geomechfile << "/" << std::endl << std::endl;

  // geomechfile << "GM_EFRAC_COHESION" << std::endl;
  // for (const auto & efrac : data.vEfrac)
  //   geomechfile << efrac.cohesion << std::endl;
  // geomechfile << "/" << std::endl << std::endl;

  // geomechfile << "GM_EFRAC_FRICTION" << std::endl;
  // for (const auto & efrac : data.vEfrac)
  //   geomechfile << efrac.friction_angle << std::endl;
  // geomechfile << "/" << std::endl << std::endl;

  // geomechfile << "GM_EFRAC_DILATION" << std::endl;
  // for (const auto & efrac : data.vEfrac)
  //   geomechfile << efrac.dilation_angle << std::endl;
  // geomechfile << "/" << std::endl << std::endl;

  // geomechfile.close();
}


void OutputDataGPRS::save_geomechanics_boundary_conditions_() const
{
  const std::string file_name = _output_path + "/" + _config.bcond_file;
  std::cout << "saving " << file_name << std::endl;

  std::ofstream out(file_name);
  // save dirichlet BC's
  if ( !_data.dirichlet_indices[0].empty() )
    save_dirichlet_component_vertices(0, "X", out);
  if ( !_data.dirichlet_indices[1].empty() )
    save_dirichlet_component_vertices(1, "Y", out);
  if ( !_data.dirichlet_indices[2].empty() )
    save_dirichlet_component_vertices(2, "Z", out);

  if ( !_data.neumann_face_indices.empty() )
  {
    out << "GMFACE_TRACTION_TXYZ" << "\n";

    for (std::size_t i=0; i<_data.neumann_face_indices.size(); ++i)
    {
      out << _data.mech_cell_numbering->face_dof(_data.neumann_face_indices[i]) + 1 << "\t";
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

void OutputDataGPRS::saveDiscreteFractureProperties(const std::string file_name)
{
  // std::cout << "write discrete fracs" << std::endl;

  // std::ofstream geomechfile;
  // geomechfile.open(file_name.c_str());
  // set<int>::iterator itsetint;

  // cout << "write all fractured faces\n";
  // geomechfile << "GMFACE_FRACTURE_TO_FLOWCELL" << std::endl;
  // for (const auto & face_it : data.dfm_faces)
  // {
  //   // geomechfile << face_it.second.ifracture + 1 << "\t";
  //   geomechfile << face_it.second.nface + 1 << "\t";
  //   if (face_it.second.coupled)
  //     geomechfile << face_it.second.nfluid + 1 << std::endl;
  //   else
  //     geomechfile << -1 << std::endl;
  // }
  // geomechfile << "/" << std::endl << std::endl;

  // // geomechfile << "GMFACE_FRACTURE_TO_FLOWCELL\n";
  // // for (const auto facet_it : data.dfm_faces)
  // // {
  // //   geomechfile << facet_it.second.nface + 1 << "\t";
  // //   geomechfile << facet_it.second.nfluid + 1 << endl;
  // // }
  // // for(itsetint = pSim->setIdenticalInternalMarker.begin();
  // //     itsetint != pSim->setIdenticalInternalMarker.end();
  // //     itsetint++, nFractures_++)
  // // {
  // //   for (int i = 0; i < pSim->nPhysicalFacets; i++)
  // //   {
  // //     if( pSim->vsPhysicalFacet[i].nmark == *itsetint )
  // //     {
  // //       geomechfile << pSim->vsPhysicalFacet[i].nface + 1 << "\t";
  // //       geomechfile << pSim->vsPhysicalFacet[i].nfluid + 1 << endl;
  // //       if( pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].nNeighbors !=2 )
  // //       {
  // //         cout << "Fracture interface # " << nFractures_ << endl;
  // //         cout << "Global interface   # " << pSim->vsPhysicalFacet[i].nface << endl;
  // //         cout << "Number od neighbors  " << pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].nNeighbors << endl;
  // //         cout << "Wrong msh file. Mesh verticies are not connected on fracture interface" << endl;
  // //         exit(0);
  // //       }
  // //     }
  // //   }
  // // }
  // // geomechfile << "/" << std::endl << std::endl;

  // geomechfile << "GMFACE_FRACTURE_CONDUCTIVITY" << std::endl;
  // for (const auto facet_it : data.dfm_faces)
  //   geomechfile << facet_it.second.conductivity << std::endl;
  // geomechfile << "/" << std::endl << std::endl;

  // geomechfile << "GMFACE_FRACTURE_REGION" << std::endl;
  // for (const auto facet_it : data.dfm_faces)
  //   geomechfile << 1 << std::endl;
  // geomechfile << "/" << std::endl << std::endl;

  // geomechfile << "GMFACE_FRACTURE_GROUP" << std::endl;
  // for (const auto facet_it : data.dfm_faces)
  //   geomechfile << 1 << std::endl;
  // geomechfile << "/" << std::endl << std::endl;

  // geomechfile.close();
}


void OutputDataGPRS::saveWells(const std::string file_name) const
{
  std::ofstream out;
  out.open(file_name.c_str());

  out << "WELSPECS" << std::endl;
  for (const auto & well : _data.wells)
  {
    out << well.name << "\tGROUP1\t";
    // connected volume + j + k empty
    assert( !well.connected_volumes.empty() );
    out << well.connected_volumes[0] << " 1 1 ";
    // reference depth
    out << -well.reference_depth << " /" << std::endl;
  }
  out << "/" << std::endl << std::endl;

  out << "COMPDAT" << std::endl;
  out << "-- name\tcell\tidk\tidk\tidk\topen\tidk\tWI\trad" << std::endl;
  for (const auto & well : _data.wells)
  {
    for (std::size_t i=0; i<well.connected_volumes.size(); ++i)
    {
      out << well.name << "\t";
      assert( !well.connected_volumes.empty() );
      out << well.connected_volumes[i] + 1 << "\t";
      // j, k1:k2 open sat_table_number
      out << "1\t1\t1\tOPEN\t1*\t";
      out << well.indices[i] * transmissibility_conversion_factor << "\t";
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
  // std::cout << "save geometry" << std::endl;
  save_geometry_();
  // save computed element data: shape functions and gradients, gauss weights, JxW
  save_fem_data_();

  save_geomechanics_keywords_();

  save_geomechanics_boundary_conditions_();
}

void OutputDataGPRS::save_cell_geometry_(std::ofstream & out, const mesh::Mesh & grid) const
{
  out << "GMCELL_NODES" << "\n";
  for (auto cell=grid.begin_active_cells(); cell!=grid.end_active_cells(); ++cell)
  {
    const auto & vertices = cell->vertices();
    out << vertices.size() << "\t";

    switch (cell->vtk_id())
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

  out << "/" << "\n\n";

  // Cell types //
  out << "GMCELL_TYPE" << "\n";
  for (auto cell=grid.begin_active_cells(); cell!=grid.end_active_cells(); ++cell)
    out <<  cell->vtk_id() << "\n";
  out << "/" << "\n\n";
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
  std::cout << "saving " << file_name << std::endl;
  std::ofstream out;
  out.open(file_name.c_str());
  // save cell data
  const std::vector<discretization::FiniteElementData> & cells = _data.fe_cell_data;
  out << "GMCELL_GAUSS_WEIGHTS" << "\n";
  for (const auto & cell : cells)
    if (!cell.points.empty())
  {
    out << _data.mech_cell_numbering->cell_dof(cell.element_index) + 1 << "\t";
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
      for (const auto &point : cell.points)
      {
        for (const angem::Point<3,double> & grad : point.grads)
          out << grad << "\t";
        out << "\n";
      }
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
  // face gauss weights
  const auto & faces = _data.fe_face_data;  // fe values and gradients for grid faces
  if (faces.empty())
    return;

  out << "GMFACE_GAUSS_WEIGHTS" << "\n";
  for (const auto & face : faces)
    if (!face.points.empty())
    {
      out << _data.mech_cell_numbering->face_dof(face.element_index) + 1 << "\t";
      out << face.points.size() << "\t";
      for (const auto & point : face.points)
        out << point.weight << "\t";
      out << "\n";
    }
  out << "/\n\n";

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
        for (const double value : point.values)
          out << value << "\t";
      out << "\n";
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
    {
      for (const auto &point : face.points)
        for (const angem::Point<3,double> & grad : point.grads)
          out << grad << "\t";
      out << "\n";
    }
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

  out.close();
}

}
