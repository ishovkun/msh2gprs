#include "OutputDataGPRS.hpp"
#include "VTKWriter.hpp"

#include <sys/stat.h>

namespace gprs_data
{

const int n_entries_per_line = 10;
static constexpr double transmissibility_conversion_factor = 0.0085267146719160104986876640419948;


OutputDataGPRS::OutputDataGPRS(const SimData & data, const GPRSOutputConfig config)
    :
    m_data(data),
    m_config(config)
{}


void OutputDataGPRS::write_output(const std::string & output_path) const
{
  // save flow data
  save_flow_data_(output_path + "/" + m_config.flow_cv_file,
                  output_path + "/" + m_config.flow_connection_file);
  // std::cout << "save geometry" << std::endl;
  // saveGeometry(output_path);

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

  if (!m_data.wells.empty())
  {
    std::cout << "save wells" << std::endl;
    saveWells(output_path + "/" + m_config.wells_file);
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

    const auto & cvs = m_data.cv_data;
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
    for (std::size_t ivar=0; ivar<m_data.output_flow_properties.size(); ++ivar)
    {
      const size_t prop_key = m_data.output_flow_properties[ivar];
      const std::string keyword = m_data.property_names[prop_key];
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

    const auto & cons = m_data.flow_connection_data;

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

void OutputDataGPRS::saveGeometry(const std::string & output_path)
{
  // // vtk output
  // // const std::string vtk_file = output_path + "reservoir_mesh.vtk";
  // // std::cout << "Saving reservoir mesh file: " << vtk_file << std::endl;
  // // IO::VTKWriter::write_vtk(grid.vertices.points, grid.cells,
  // //                          grid.shape_ids, vtk_file);

  // // if (data.dfm_faces.size() > 0)
  // // { // DFM frac geometry
  //   // const std::string vtk_dfm_file = output_path + "dfm.vtk";
  //   // IO::VTKWriter::write_geometry(data.dfm_master_grid.get_vertices(),
  //   //                               data.dfm_master_grid.get_polygons(),
  //   //                               vtk_dfm_file);

  // // }

  // // gprs output
  // std::string outstring;
  // std::ofstream geomechfile;

  // // GEOMETRY
  // outstring = output_path + data.config.mechanics_domain_file;
  // std::cout << "writing file " << outstring << std::endl;

  // geomechfile.open(outstring.c_str());
  // geomechfile << "GMDIMS" << endl;

  // geomechfile << grid.n_vertices() << "\t"
  //             << grid.n_cells() << "\t"
  //             << grid.n_faces();
  // geomechfile << "/" << endl << endl;

  // geomechfile.precision(6);
  // cout << "write all coordinates\n";
  // geomechfile << "GMNODE_COORDS" << endl;
  // for (const auto & vertex : grid.get_vertices())
  //     geomechfile << vertex[0] << "\t"
  //                 << vertex[1] << "\t"
  //                 << vertex[2] << "\n";
  // geomechfile << "/" << std::endl << std::endl;

  // cout << "write all elements\n";
  // geomechfile << "GMCELL_NODES" << endl;
  // for (auto cell=grid.begin_cells(); cell!=grid.end_cells(); ++cell)
  // {
  //   const auto & vertices = cell.vertex_indices();
  //   geomechfile << vertices.size() << "\t";

  //   switch (cell.vtk_id())
  //   {
  //     case 25: // super wierd element 25
  //       {
  //         for (int j = 0; j < 8; j++)
  //           geomechfile << vertices[j] + 1 << "\t";
  //         geomechfile << vertices[8] + 1 << "\t";
  //         geomechfile << vertices[11] + 1 << "\t";
  //         geomechfile << vertices[13] + 1 << "\t";
  //         geomechfile << vertices[9] + 1 << "\t";

  //         geomechfile << vertices[16] + 1 << "\t";
  //         geomechfile << vertices[18] + 1 << "\t";
  //         geomechfile << vertices[19] + 1 << "\t";
  //         geomechfile << vertices[17] + 1 << "\t";

  //         geomechfile << vertices[10] + 1 << "\t";
  //         geomechfile << vertices[12] + 1 << "\t";
  //         geomechfile << vertices[14] + 1 << "\t";
  //         geomechfile << vertices[15] + 1 << "\t";
  //         break;
  //       }
  //     case 26:
  //       {
  //         for (int j = 0; j < 6; j++)
  //           geomechfile << vertices[j] + 1 << "\t";

  //         geomechfile << vertices[6] + 1 << "\t";
  //         geomechfile << vertices[9] + 1 << "\t";
  //         geomechfile << vertices[7] + 1 << "\t";

  //         geomechfile << vertices[12] + 1 << "\t";
  //         geomechfile << vertices[14] + 1 << "\t";
  //         geomechfile << vertices[13] + 1 << "\t";

  //         geomechfile << vertices[8] + 1 << "\t";
  //         geomechfile << vertices[10] + 1 << "\t";
  //         geomechfile << vertices[11] + 1 << "\t";
  //         break;
  //       }
  //     default:
  //       {
  //         for (const auto vertex : vertices)
  //           geomechfile << vertex + 1 << "\t";
  //         break;
  //       }
  //   }

  //   geomechfile << std::endl;
  // }

  // geomechfile << "/" << endl << endl;

  // geomechfile << "GMCELL_TYPE" << endl;
  // for (auto cell=grid.begin_cells(); cell!=grid.end_cells(); ++cell)
  //   geomechfile <<  cell.vtk_id() << std::endl;
  // geomechfile << "/" << std::endl << std::endl;

  // geomechfile << "GMCELL_TO_FLOWCELLS" << endl;
  // for (const auto & flow_cells : data.gm_cell_to_flow_cell)
  // {
  //   const std::size_t n_connected_elements = flow_cells.size();
  //   if (n_connected_elements == 0)
  //     geomechfile << 1 << "\t" << -1 << std::endl;
  //   else
  //   {
  //     geomechfile << n_connected_elements << "\t";
  //     for (const std::size_t ielement : flow_cells)
  //       geomechfile << ielement + 1 << "\t";
  //     geomechfile << std::endl;
  //   }
  // }
  // geomechfile << "/\n\n";

  // std::cout << "write all faces\n";
  // geomechfile << "GMFACE_NODES\n";
  // for (auto face=grid.begin_faces(); face!=grid.end_faces(); ++face)
  // {
  //   const std::vector<std::size_t> ivertices = face.vertex_indices();
  //   geomechfile << ivertices.size() << "\t";
  //   for (const std::size_t ivertex : ivertices)
  //     geomechfile << ivertex + 1 << "\t";
  //   geomechfile << std::endl;
  // }

  // geomechfile << "/" << std::endl << std::endl;

  // geomechfile << "GMFACE_TYPE" << std::endl;
  // for (auto face=grid.begin_faces(); face!=grid.end_faces(); ++face)
  //   geomechfile << face.vtk_id() << std::endl;
  // geomechfile << "/" << std::endl << std::endl;

  // std::cout << "writing face-cell connection" << std::endl;
  // geomechfile << "GMFACE_GMCELLS" << std::endl;
  // for (auto face=grid.begin_faces(); face!=grid.end_faces(); ++face)
  // {
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

  //       geomechfile << neighbors.size() << "\t";
  //       for (const auto & neighbor : neighbors)
  //         geomechfile << neighbor + 1 << "\t";
  //       geomechfile << std::endl;
  //     }
  //   }
  //   else
  //   {
  //     const auto & neighbors = face.neighbors();
  //     geomechfile << neighbors.size() << "\t";
  //     for (const auto & neighbor : neighbors)
  //       geomechfile << neighbor + 1 << "\t";
  //     geomechfile << std::endl;
  //   }
  // }
  // geomechfile << "/" << std::endl << std::endl;

}


void OutputDataGPRS::saveGeomechDataNewKeywords(const std::string file_name)
{
  // std::cout << "Writing all domain from user input properties" << std::endl;
  // {  // write domain properties
  //   std::ofstream geomechfile;
  //   geomechfile.open(file_name.c_str());

  //   for (std::size_t ivar=0; ivar<data.rockPropNames.size(); ++ivar)
  //   {
  //     if ( data.config.expression_type[ivar] != 1 )  // only mechanics kwds
  //       continue;

  //     geomechfile << data.rockPropNames[ivar] << std::endl;
  //     for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
  //     {
  //       const std::size_t icell = cell->index();
  //       geomechfile << data.vsCellRockProps[icell].v_props[ivar] << "\t";
  //       if ((icell + 1) % n_entries_per_line == 0)
  //         geomechfile << std::endl;
  //     }
  //     geomechfile << std::endl << "/" << std::endl << std::endl;
  //   }

  //   geomechfile.close();
  // }
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


void OutputDataGPRS::saveBoundaryConditions(const std::string file_name)
{
  // std::ofstream geomechfile;
  // std::cout << "Computing Dirichlet nodes" << std::endl;

  // /* this chunk of code find dirichlet nodes
  //  * Algorithm:
  //  * 1. First add all Dirichlet nodes specified by the user
  //  * to the set
  //  * 2. loop Dirichlet faces.
  //  * loop nodes withing dirichlet faces
  //  * add node as a part of dirichlet face
  //  */
  // const int dim = 3;
  // std::vector<std::unordered_set<std::size_t>> setDisp(dim);
  // std::vector<std::vector<std::size_t>> vv_disp_node_ind(dim);
  // std::vector<std::vector<double>> vv_disp_node_values(dim);

  // // search tolerance
  // const double tol = data.config.node_search_tolerance;

  // for (std::size_t ivert=0; ivert<grid.n_vertices(); ++ivert)
  // {
  //   const auto & vertex = grid.vertex(ivert);

  //   for (const auto & bc_node : data.config.bc_nodes)
  //   {
  //     if (bc_node.coord.distance(vertex) < tol)
  //     {
  //       for (int d=0; d<dim; ++d)
  //         if (bc_node.value[d] != data.config.nan)
  //         {
  //           setDisp[d].insert(ivert);
  //           vv_disp_node_values[d].push_back(bc_node.value[d]);
  //           vv_disp_node_ind[d].push_back(ivert);
  //         }
  //     }
  //   }
  // }

  // if ( data.n_dirichlet_faces > 0 )
  // for (auto face=grid.begin_faces(); face!=grid.end_faces(); ++face)
  //     if (data.is_boundary(face->marker()))
  //     {
  //       const auto facet_it = data.boundary_faces.find(face->master_index());
  //       if (facet_it != data.boundary_faces.end())
  //       {
  //         if (facet_it->second.ntype == 1) // dirichlet
  //           for (std::size_t d=0; d<dim; ++d)
  //             if (facet_it->second.condition[d] != data.config.nan)
  //               for (const auto ivertex : face->vertices())
  //               {
  //                 const auto ret = setDisp[d].insert(ivertex);
  //                 // check if already in set
  //                 if(ret.second)
  //                 {
  //                   vv_disp_node_ind[d].push_back(ivertex);
  //                   vv_disp_node_values[d].push_back(facet_it->second.condition[d]);
  //                 } // end if vertex is new
  //               }
  //       }
  //       else  // should not happen
  //       {
  //         std::cout << "face " << face->index()
  //                   << " is not a boundary face!! Aborting!!!"
  //                   << std::endl;
  //         abort();
  //       }
  //     }

  // setDisp.clear();

  // geomechfile.open(file_name.c_str());

  // // save dirichlet faces
  // for (std::size_t j=0; j<dim; ++j)
  //   if (vv_disp_node_ind[j].size() > 0)
  //   {
  //     std::cout << "write all Dirichlet nodes" << std::endl;
  //     switch (j)
  //     {
  //       case 0:
  //         geomechfile << "GMNODE_BCDISPX" << std::endl;
  //         break;
  //       case 1:
  //         geomechfile << "GMNODE_BCDISPY" << std::endl;
  //         break;
  //       case 2:
  //         geomechfile << "GMNODE_BCDISPZ" << std::endl;
  //         break;
  //     }

  //     for(std::size_t i = 0; i < vv_disp_node_ind[j].size(); ++i)
  //     {
  //       geomechfile << vv_disp_node_ind[j][i] + 1 << "\t";
  //       geomechfile << vv_disp_node_values[j][i] << std::endl;
  //     }
  //     geomechfile << "/" << std::endl << std::endl;
  //   }

  // if ( data.n_neumann_faces > 0 )
  // {
  //   std::cout << "write all Neumann faces" << std::endl;

  //   geomechfile << "GMFACE_TRACTION_TXYZ" << std::endl;
  //   for (const auto & facet_it : data.boundary_faces)
  //     if (facet_it.second.ntype == 2)  // neumann
  //     {
  //       geomechfile << facet_it.second.nface + 1 << "\t";
  //       geomechfile << facet_it.second.condition << std::endl;
  //     }

  //   geomechfile << "/" << std::endl << std::endl;
  // }

  // geomechfile.close();
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
  for (const auto & well : m_data.wells)
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
  for (const auto & well : m_data.wells)
  {
    for (std::size_t i=0; i<well.connected_volumes.size(); ++i)
    {
      out << well.name << "\t";
      assert( !well.connected_volumes.empty() );
      out << well.connected_volumes[i] + 1 << "\t";
      // j, k1:k2 open sat_table_number
      out << "1\t1\t1\tOPEN\t1*\t";
      out << well.indices[i] * transmissibility_conversion_factor << "\t";
      out << 2*well.radius << "\t";
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

}
