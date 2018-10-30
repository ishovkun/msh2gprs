#include "femout.hpp"
#include "VTKWriter.hpp"

#include <sys/stat.h>

const int n_entries_per_line = 10;


OutputData::OutputData(SimData * pSimData)
{
  pSim = pSimData;
}

OutputData::~OutputData()
{
}


void OutputData::write_output(const std::string & output_path)
{
  writeGeomechDataNewKeywords(output_path);
  CalcTranses::save_output(pSim->flow_data, output_path);
}


void OutputData::writeGeomechDataNewKeywords(const std::string & output_path)
{
  stringstream out;
  out << "model";
  out << "/";
  // string outputPath_ = ""; // out.str();

// #if 0 //defined(_WIN32)
//   _mkdir(output_path.c_str()); // create directory
// //#else
//   mkdir(outputPath_.c_str(), 0777);
// #endif
  // string outstring = pSim->outstream;
  string outstring;

  ofstream geomechfile;

  outstring = output_path + "gm_depth.txt";
  std::cout << "writing file: " << outstring  << std::endl;
  geomechfile.open(outstring.c_str());
  geomechfile << "GDEPTH" << endl;

  for (int i = 0; i < pSim->nCells; i++)
  {
	  geomechfile << pSim->vsCellCustom[i].center [2] << "\n";
  }
  geomechfile << "/" << endl << endl;
  geomechfile.close();

  // GEOMETRY
  outstring =   output_path + "gm_geometry.txt";
  geomechfile.open(outstring.c_str());
  geomechfile << "GMDIMS" << endl;
  geomechfile <<  pSim->nNodes << "\t" << pSim->nCells << "\t" << pSim->nFaces;
  geomechfile << "/" << endl << endl;

  const std::string vtk_file = output_path + "reservoir_mesh.vtk";
  IO::VTKWriter::write_vtk(pSim->vvVrtxCoords,
                           pSim->vsCellCustom,
                           vtk_file);

  geomechfile.precision(12);
  cout << "write all coordinates\n";
  geomechfile << "GMNODE_COORDS" << endl;
  for (int i = 0 ; i < pSim->nNodes; i++)
  {
    geomechfile << pSim->vvVrtxCoords[i][0] << "\t"<< pSim->vvVrtxCoords[i][1] << "\t"<< pSim->vvVrtxCoords[i][2] << "\n";
  }
  geomechfile << "/" << endl << endl;

  cout << "write all elements\n";
  geomechfile << "GMCELL_NODES" << endl;
  for (int i = 0; i < pSim->nCells; i++)
  {
    geomechfile << pSim->vsCellCustom[i].nVertices << "\t";
    if(pSim->vsCellCustom[i].vtkIndex == 25)
    {
      // super wierd element 25
      for (int j = 0; j < 8; j++)
        geomechfile << pSim->vsCellCustom[i].vVertices[j] + 1 << "\t";

      geomechfile << pSim->vsCellCustom[i].vVertices[8] + 1 << "\t";
      geomechfile << pSim->vsCellCustom[i].vVertices[11] + 1 << "\t";
      geomechfile << pSim->vsCellCustom[i].vVertices[13] + 1 << "\t";
      geomechfile << pSim->vsCellCustom[i].vVertices[9] + 1 << "\t";

      geomechfile << pSim->vsCellCustom[i].vVertices[16] + 1 << "\t";
      geomechfile << pSim->vsCellCustom[i].vVertices[18] + 1 << "\t";
      geomechfile << pSim->vsCellCustom[i].vVertices[19] + 1 << "\t";
      geomechfile << pSim->vsCellCustom[i].vVertices[17] + 1 << "\t";

      geomechfile << pSim->vsCellCustom[i].vVertices[10] + 1 << "\t";
      geomechfile << pSim->vsCellCustom[i].vVertices[12] + 1 << "\t";
      geomechfile << pSim->vsCellCustom[i].vVertices[14] + 1 << "\t";
      geomechfile << pSim->vsCellCustom[i].vVertices[15] + 1 << "\t";
      geomechfile << endl;
    }
    else if(pSim->vsCellCustom[i].vtkIndex == 26)
    {
      for (int j = 0; j < 6; j++)
        geomechfile << pSim->vsCellCustom[i].vVertices[j] + 1 << "\t";

      geomechfile << pSim->vsCellCustom[i].vVertices[6] + 1 << "\t";
      geomechfile << pSim->vsCellCustom[i].vVertices[9] + 1 << "\t";
      geomechfile << pSim->vsCellCustom[i].vVertices[7] + 1 << "\t";

      geomechfile << pSim->vsCellCustom[i].vVertices[12] + 1 << "\t";
      geomechfile << pSim->vsCellCustom[i].vVertices[14] + 1 << "\t";
      geomechfile << pSim->vsCellCustom[i].vVertices[13] + 1 << "\t";

      geomechfile << pSim->vsCellCustom[i].vVertices[8] + 1 << "\t";
      geomechfile << pSim->vsCellCustom[i].vVertices[10] + 1 << "\t";
      geomechfile << pSim->vsCellCustom[i].vVertices[11] + 1 << "\t";
      geomechfile << endl;
    }
    else
    {
      for (int j = 0; j < pSim->vsCellCustom[i].nVertices; j++)
        geomechfile << pSim->vsCellCustom[i].vVertices[j] + 1 << "\t";
      geomechfile << endl;
    }
  }
  geomechfile << "/" << endl << endl;


  geomechfile << "GMCELL_TYPE" << endl;
  for (int i = 0; i < pSim->nCells; i++)
    geomechfile <<  pSim->vsCellCustom[i].vtkIndex << endl;
  geomechfile << "/" << endl << endl;

  geomechfile << "GMCELL_TO_FLOWCELLS" << endl;
  int i_ = 0;
  for (int i = 0; i < pSim->nCells; i++)
  {
    // geomechfile << 1 << "\t" << pSim->vsCellCustom[i].fluidElement + 1 << endl;
    geomechfile << 1 << "\t" << -1 << endl;
  }
  geomechfile << "/\n\n";

  cout << "write all faces\n";
  geomechfile << "GMFACE_NODES\n";
  for (int i = 0; i < pSim->nFaces; i++)
  {
    geomechfile << pSim->vsFaceCustom[i].nVertices << "\t";

    for (int j = 0; j < pSim->vsFaceCustom[i].nVertices; j++) geomechfile << pSim->vsFaceCustom[i].vVertices[j] + 1 << "\t";

    geomechfile << endl;
  }
  geomechfile << "/" << endl << endl;

  geomechfile << "GMFACE_TYPE" << endl;
  for (int i = 0; i < pSim->nFaces; i++)
    geomechfile <<  pSim->vsFaceCustom[i].vtkIndex << endl;
  geomechfile << "/" << endl << endl;

  geomechfile << "GMFACE_GMCELLS" << endl;
  for (int i = 0; i < pSim->nFaces; i++)
  {
      geomechfile << pSim->vsFaceCustom[i].nNeighbors << "\t";

      for (int j = 0; j < pSim->vsFaceCustom[i].nNeighbors; j++)
        geomechfile << pSim->vsFaceCustom[i].vNeighbors[j] + 1 << "\t";

      geomechfile << endl;
  }
  geomechfile << "/" << endl << endl;
  geomechfile.close();

  std::cout << "Writing all domain from user input properties" << std::endl;
  {  // write domain properties
    outstring =   output_path + pSim->config.domain_file;
    geomechfile.open(outstring.c_str());

    for (std::size_t ivar=0; ivar<pSim->rockPropNames.size(); ++ivar)
    {
      if ( pSim->config.expression_type[ivar] != 1 )  // only mechanics kwds
        continue;

      geomechfile << pSim->rockPropNames[ivar] << std::endl;
      for (std::size_t icell=0; icell<pSim->nCells; ++icell)
      {
        geomechfile << pSim->vsCellRockProps[icell].v_props[ivar] << "\t";
        if ((icell + 1) % n_entries_per_line == 0)
          geomechfile << std::endl;
      }
      geomechfile << std::endl << "/" << std::endl << std::endl;
    }

    geomechfile.close();
  }

  // Strong discontinuity
  if (pSim->vEfrac.size() > 0)
  {
    std::cout  << "Writing SDA props" << std::endl;
    outstring =   output_path + "gm_SDA.txt";
    geomechfile.open(outstring.c_str());

    geomechfile << "GM_EFRAC_CELLS" << std::endl;
    for (const auto & efrac : pSim->vEfrac)
    {
      geomechfile << efrac.cells.size() << std::endl << "\t";
      for (std::size_t i=0; i<efrac.cells.size(); ++i)
      {
        geomechfile << efrac.cells[i] + 1 << "\t";
        if ((i+1) % n_entries_per_line == 0)
          geomechfile << std::endl;
        if (i == efrac.cells.size() - 1)
          geomechfile << std::endl;
      }
    }
    geomechfile << "/" << std::endl << std::endl;

    geomechfile << "GM_EFRAC_POINTS" << std::endl;
    for (const auto & efrac : pSim->vEfrac)
      for (std::size_t i=0; i<efrac.points.size(); ++i)
          geomechfile << efrac.points[i] << std::endl;
    geomechfile << "/" << std::endl << std::endl;

    geomechfile << "GM_EFRAC_DIP" << std::endl;
    for (const auto & efrac : pSim->vEfrac)
      for (std::size_t i=0; i<efrac.points.size(); ++i)
      {
        geomechfile << efrac.dip[i] << "\t";
        if ((i+1) % n_entries_per_line == 0)
          geomechfile << std::endl;
      }
    geomechfile << "/" << std::endl << std::endl;

    geomechfile << "GM_EFRAC_STRIKE" << std::endl;
    for (const auto & efrac : pSim->vEfrac)
      for (std::size_t i=0; i<efrac.points.size(); ++i)
      {
        geomechfile << efrac.strike[i] << "\t";
        if ((i+1)%n_entries_per_line == 0) geomechfile << std::endl;
      }
    geomechfile << "/" << std::endl << std::endl;

    geomechfile << "GM_EFRAC_COHESION" << std::endl;
    for (const auto & efrac : pSim->vEfrac)
      geomechfile << efrac.cohesion << std::endl;
    geomechfile << "/" << std::endl << std::endl;

    geomechfile << "GM_EFRAC_FRICTION" << std::endl;
    for (const auto & efrac : pSim->vEfrac)
      geomechfile << efrac.friction_angle << std::endl;
    geomechfile << "/" << std::endl << std::endl;

    geomechfile << "GM_EFRAC_DILATION" << std::endl;
    for (const auto & efrac : pSim->vEfrac)
      geomechfile << efrac.dilation_angle << std::endl;
    geomechfile << "/" << std::endl << std::endl;

    geomechfile.close();

    // write vtk data
    // std::cout  << "Writing SDA props" << std::endl;
    {
      std::size_t n_efrac_vertices = 0;
      // make up a vector of all sda vertices
      for (const auto & efrac : pSim->vEfrac)
        n_efrac_vertices += efrac.mesh.vertices.size();

      std::vector<angem::Point<3,double>> efrac_verts(n_efrac_vertices);
      std::size_t ivertex = 0;
      for (const auto & efrac : pSim->vEfrac)
        for (const auto & p : efrac.mesh.vertices)
        {
          efrac_verts[ivertex] = p;
          ivertex++;
        }

      std::size_t n_efrac_elements = 0;
      std::size_t vind_size_total = 0;

      for (const auto & efrac : pSim->vEfrac)
      {
        // n_efrac_elements += efrac.vIndices.size();
        // for (const auto & vec : efrac.vIndices)
        //   vind_size_total += vec.size();
        n_efrac_elements += efrac.mesh.polygons.size();
        for (const auto & poly : efrac.mesh.polygons)
          vind_size_total += poly.size();
      }

      std::vector<std::vector<std::size_t>> efrac_cells(n_efrac_elements);
      std::size_t ielement = 0;
      std::size_t shift = 0;
      for (const auto & efrac : pSim->vEfrac)
      {
        for (const auto & cell : efrac.mesh.polygons)
        {
          efrac_cells[ielement].resize(cell.size());
          for (short v=0; v<cell.size(); ++v)
            efrac_cells[ielement][v] = shift + cell[v];

          ielement++;
        }
        // shift += efrac.vVertices.size();
        shift += efrac.mesh.vertices.size();
      }

      const std::string vtk_file = output_path + "efrac.vtk";
      IO::VTKWriter::write_vtk(efrac_verts, efrac_cells, vtk_file);

      // const std::string vtk_file2 = output_path + "ababa.vtk";
      // IO::VTKWriter::write_vtk(pSim->vEfrac[0].mesh.vertices.points,
      //                          pSim->vEfrac[0].mesh.polygons,
      //                          vtk_file2);
    }
    // outstring =   output_path + "efrac.vtk";
    // geomechfile.open(outstring.c_str());

    // // const auto & efrac = pSim->vEfrac[0];
    // geomechfile << "# vtk DataFile Version 2.0 \n";
    // geomechfile << "3D Fractures \n";
    // geomechfile << "ASCII \n \n";
    // geomechfile << "DATASET UNSTRUCTURED_GRID \n";

    // std::size_t n_points = 0;
    // for (const auto & efrac : pSim->vEfrac)
    //   n_points += efrac.vVertices.size();
    // geomechfile << "POINTS" << "\t"
    //             << n_points << " float"
    //             << std::endl;

    // for (const auto & efrac : pSim->vEfrac)
    //   for (const auto & p : efrac.vVertices)
    //     geomechfile << p << std::endl;
    // geomechfile << std::endl;

    // // count number of entries in vindices
    // std::size_t vind_size_total = 0;
    // std::size_t n_cells = 0;
    // for (const auto & efrac : pSim->vEfrac)
    // {
    //   n_cells += efrac.vIndices.size();
    //   for (const auto & vec : efrac.vIndices)
    //     vind_size_total += vec.size();
    // }

    // geomechfile << "CELLS" << "\t"
    //             << n_cells << "\t"
    //             << vind_size_total + n_cells
    //             << std::endl;

    // std::size_t shift = 0;
    // for (const auto & efrac : pSim->vEfrac)
    // {
    //   for (const auto & cell : efrac.vIndices)
    //   {
    //     geomechfile << cell.size() << "\t";
    //     for (const std::size_t & ivert : cell)
    //       geomechfile << shift + ivert << "\t";
    //     geomechfile << std::endl;

    //   }
    //   shift += efrac.vVertices.size();
    // }

    // geomechfile << std::endl;
    // geomechfile << "CELL_TYPES" << "\t" << n_cells << std::endl;
    // for (const auto & efrac : pSim->vEfrac)
    //   for (const auto & cell : efrac.vIndices)
    //   {
    //     // vtk indices
    //     if (cell.size() == 4)  // quad
    //       geomechfile << 9 << std::endl;
    //     else if (cell.size() == 3)  // triangle
    //       geomechfile << 5 << std::endl;
    //     else if (cell.size() == 5 || cell.size() == 6)  // polygon
    //       geomechfile << 7 << std::endl;
    //     else
    //     {
    //       std::cout << "unknown cell type : " << cell.size() << " vertices" << std::endl;
    //       exit(-1);
    //     }
    //   }

    // geomechfile.close();
  }  // end if efracs exist


  outstring =   output_path + "gm_bcond.txt";
  geomechfile.open(outstring.c_str());

  if ( pSim->nNeumannFaces > 0 )
  {
    cout << "write all Neumann faces\n";

    geomechfile << "GMFACE_TRACTION_TXYZ" << endl;
    for ( std::size_t i = 0; i < pSim->nPhysicalFacets; i++ )
      if (pSim->vsPhysicalFacet[i].ntype == 2)
      {
        geomechfile << pSim->vsPhysicalFacet[i].nface + 1 << "\t";
        geomechfile << pSim->vsPhysicalFacet[i].condition << endl;
      }
    geomechfile << "/\n\n";
  }


  std::cout << "Computing Dirichlet nodes" << std::endl;
  /* this chunk of code find dirichlet nodes
   * Algorithm:
   * 1. First add all Dirichlet nodes specified by the user
   * to the set
   * 2. loop Dirichlet faces.
   * loop nodes withing dirichlet faces
   * add node as a part of dirichlet face
   */
  const int dim = 3;
  std::vector<set<int>> setDisp(dim);
  std::vector<std::vector<std::size_t>> vv_disp_node_ind(dim);
  std::vector<std::vector<double>> vv_disp_node_values(dim);

  // search tolerance
  const double tol = pSim->config.node_search_tolerance;

  for (std::size_t ivert=0; ivert<pSim->vvVrtxCoords.size(); ++ivert)
  {
    const auto & vertex = pSim->vvVrtxCoords[ivert];

    for (const auto & bc_node : pSim->config.bc_nodes)
    {
      if (bc_node.coord.distance(vertex) < tol)
      {
        for (int d=0; d<dim; ++d)
          if (bc_node.value[d] != pSim->config.nan)
          {
            setDisp[d].insert(ivert);
            vv_disp_node_values[d].push_back(bc_node.value[d]);
            vv_disp_node_ind[d].push_back(ivert);
          }
      }
    }
  }

  if ( pSim->nDirichletFaces > 0 )
    for ( int i = 0; i < pSim->nPhysicalFacets; i++ )
      if ( pSim->vsPhysicalFacet[i].ntype == 1 )  // dirichlet
        for (std::size_t d=0; d<dim; ++d)
          if ( pSim->vsPhysicalFacet[i].condition[d] != pSim->config.nan )
          {
            const int nvrtx = pSim->vsFaceCustom
                [ pSim->vsPhysicalFacet[i].nface ].nVertices;
            for ( int ivrtx = 0; ivrtx < nvrtx; ++ivrtx )
            {
              const std::size_t global_vertex_index =
                  pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].vVertices[ivrtx];
              const auto ret = setDisp[d].insert(global_vertex_index);

              // check if already in set
              if(ret.second)
              {
                vv_disp_node_ind[d].push_back(global_vertex_index);
                vv_disp_node_values[d].push_back(pSim->vsPhysicalFacet[i].condition[d]);
              } // end if vertex is new
            }  // end face vertices loop
          } // end nan condition

  setDisp.clear();

  std::cout << "write all Dirichlet nodes" << std::endl;
  for (std::size_t j=0; j<dim; ++j)
    if (vv_disp_node_ind[j].size() > 0)
    {
      switch (j)
      {
        case 0:
          geomechfile << "GMNODE_BCDISPX" << std::endl;
          break;
        case 1:
          geomechfile << "GMNODE_BCDISPY" << std::endl;
          break;
        case 2:
          geomechfile << "GMNODE_BCDISPZ" << std::endl;
          break;
      }
      for(int i = 0; i < vv_disp_node_ind[j].size(); ++i)
      {
        geomechfile <<  vv_disp_node_ind[j][i] + 1 << "\t";
        geomechfile << vv_disp_node_values[j][i] << std::endl;
      }
      geomechfile << "/" << std::endl << std::endl;
    }

  geomechfile.close();

  // std::cout << "write fracs" << std::endl;

  // outstring =   output_path + "gm_fractures.txt";
  // geomechfile.open(outstring.c_str());
  // set<int>::iterator itsetint;

 // int counter = 0;
 // if (pSim->nInternalBoundaryFaces > 0)
 // {
 //    int nFractures_ = 0;
 //    cout << "write all fractured faces\n";

 //    geomechfile << "GMFACE_FRACTURE_TO_FLOWCELL\n";
 //    for(itsetint = pSim->setIdenticalInternalMarker.begin(); itsetint != pSim->setIdenticalInternalMarker.end(); itsetint++, nFractures_++)
 //    {
 //      for (int i = 0; i < pSim->nPhysicalFacets; i++)
 //      {
 //        if( pSim->vsPhysicalFacet[i].nmark == *itsetint )
 //        {
 //    geomechfile << pSim->vsPhysicalFacet[i].nface + 1 << "\t";
 //          geomechfile << pSim->vsPhysicalFacet[i].nfluid + 1 << endl;
 //          if( pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].nNeighbors !=2 )
 //          {
 //            cout << "Fracture interface # " << nFractures_ << endl;
 //            cout << "Global interface   # " << pSim->vsPhysicalFacet[i].nface << endl;
 //            cout << "Number od neighbors  " << pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].nNeighbors << endl;
 //            cout << "Wrong msh file. Mesh verticies are not connected on fracture interface" << endl;
 //            exit(0);
 //          }
 //        }
 //      }
 //    }
 //    geomechfile << "/" << endl << endl;

 //    nFractures_ = 0;
 //    geomechfile << "GMFACE_FRACTURE_CONDUCTIVITY" << endl;
 //    for(itsetint = pSim->setIdenticalInternalMarker.begin(); itsetint != pSim->setIdenticalInternalMarker.end(); itsetint++, nFractures_++)
 //    {
 //      for (int i = 0; i < pSim->nPhysicalFacets; i++)
 //      {
 //        if( pSim->vsPhysicalFacet[i].nmark == *itsetint )
 //    geomechfile << pSim->vsFaceCustom[pSim->vsPhysicalFacet[i].nface].conductivity << endl;
 //      }
 //    }
 //    geomechfile << "/" << endl << endl;

 //    nFractures_ = 0;
 //    geomechfile << "GMFACE_FRACTURE_REGION\n";
 //    for(itsetint = pSim->setIdenticalInternalMarker.begin(); itsetint != pSim->setIdenticalInternalMarker.end(); itsetint++, nFractures_++)
 //    {
 //      for (int i = 0; i < pSim->nPhysicalFacets; i++)
 //      {
 //        if( pSim->vsPhysicalFacet[i].nmark == *itsetint )
 //        {
 //    geomechfile << 1 << endl;
 //        }
 //      }
 //    }
 //    geomechfile << "/" << endl << endl;

 //    nFractures_ = 0;
 //    geomechfile << "GMFACE_FRACTURE_GROUP\n";
 //    for(itsetint = pSim->setIdenticalInternalMarker.begin(); itsetint != pSim->setIdenticalInternalMarker.end(); itsetint++, nFractures_++)
 //    {
 //      for (int i = 0; i < pSim->nPhysicalFacets; i++)
 //      {
 //        if( pSim->vsPhysicalFacet[i].nmark == *itsetint )
 //        {
 //    geomechfile << nFractures_ + 1 << endl;
 //        }
 //      }
 //    }
 //    geomechfile << "/" << endl << endl;
 //    geomechfile.close();
 // }
 // geomechfile.close();

  // std::cout << "pressure" << std::endl;

  // outstring =   output_path + "fl_pres.txt";
  // geomechfile.open(outstring.c_str());

  // // fractures first
  // geomechfile << "PRESSURE" << endl;
  // for(int iface = 0; iface < pSim->nFaces; iface++)
  // {
  //   if( pSim->vsFaceCustom[iface].nMarker > 0)
  //     geomechfile << (pSim->vsCellRockProps[pSim->vsFaceCustom[iface].vNeighbors[0]].pressure +
  //                     pSim->vsCellRockProps[pSim->vsFaceCustom[iface].vNeighbors[1]].pressure) / 2.0 << endl;
  // }
  // // matrix
  // for(int ib = 0; ib < pSim->nCells; ++ib)
  //   geomechfile << pSim->vsCellRockProps[ib].pressure << endl;
  // geomechfile << "\n/\n\n";
  // geomechfile.close();

  // std::cout << "temperature" << std::endl;

  // outstring =   output_path + "fl_temp.txt";
  // geomechfile.open(outstring.c_str());

  // // fractures first
  // geomechfile << "RTEMP" << endl;
  // for(int iface = 0; iface < pSim->nFaces; iface++)
  // {
  //   // internal suface
  //   if( pSim->vsFaceCustom[iface].nMarker > 0)
  //   {
  //     int n1_ = pSim->vsFaceCustom[iface].vNeighbors[0];
  //     int n2_ = pSim->vsFaceCustom[iface].vNeighbors[1];
  //     geomechfile << min(pSim->vsCellRockProps[n1_].temp, pSim->vsCellRockProps[n2_].temp) << endl;
  //   }
  // }
  // // matrix
  // for(int ib = 0; ib < pSim->nCells; ++ib)
  //   geomechfile << pSim->vsCellRockProps[ib].temp << endl;

  // geomechfile << "\n/\n\n";
  // geomechfile.close();

  // std::cout << "zmf" << std::endl;

  // outstring =   output_path + "fl_zmf.txt";
  // geomechfile.open(outstring.c_str());

  // geomechfile << "ZMF\n";
  // for ( int k = 0; k < 5; ++k )
  // {
  //   for ( int iface = 0; iface < pSim->nFaces; iface++ )
  //   {
  //     // internal suface
  //     if ( pSim->vsFaceCustom[iface].nMarker > 0 )
  //     {
  //       int n1_ = pSim->vsFaceCustom[iface].vNeighbors[0];
  //       int n2_ = pSim->vsFaceCustom[iface].vNeighbors[1];
  //       geomechfile << min ( pSim->vsCellRockProps[n1_].zmf[k], pSim->vsCellRockProps[n2_].zmf[k] ) << endl;
  //     }
  //   }

  //   for(int ib = 0; ib < pSim->nCells; ++ib)
  //     geomechfile << pSim->vsCellRockProps[ib].zmf[k] << endl;
  // }

  // geomechfile << "\n/\n\n";
  // geomechfile.close();


  // cout << "write all wells\n";
  // if ( pSim->nWells > 0 )
  // {
  //   outstring =   output_path + "fl_wells.txt";
  //   geomechfile.open ( outstring.c_str() );
  //   geomechfile <<  "WELSPECS\n";

  //   for ( int iw = 0; iw < pSim->nWells; iw++ )
  //   {
  //     if(pSim->vsWell[iw].vID.size() > 0)
  //       geomechfile <<  "W" << iw
  //                   << " 1* " << pSim->vsWell[iw].vID[0]+1
  //                   << " 1 " <<  pSim->vsWell[iw].datum << " WATER /" << endl;
  //   }

  //   geomechfile << "/\n\n";

  //   geomechfile <<  "COMPDAT\n";
  //   for ( int iw = 0; iw < pSim->nWells; iw++ )
  //   {
  //     for ( int i = 0; i < pSim->vsWell[iw].vID.size(); ++i )
  //       geomechfile <<  "W" << iw << "\t" << pSim->vsWell[iw].vID[i]+1 <<" 1 1 1 OPEN 1* " << pSim->vsWell[iw].vWi[i] << " 4* Z/" << endl;
  //   }
  //   geomechfile << "/\n\n";
  //   geomechfile.close();
  // }
}
