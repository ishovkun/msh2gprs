#include "femout.hpp"
#include "VTKWriter.hpp"

#include <sys/stat.h>

const int n_entries_per_line = 10;


OutputData::OutputData(SimData & sim_data,
                       mesh::Mesh & grid)
    :
    data(sim_data),
    grid(grid)
{}

OutputData::~OutputData()
{
}


void OutputData::write_output(const std::string & output_path)
{
  std::cout << "save geometry" << std::endl;
  saveGeometry(output_path);

  std::cout << "save custom keyswords" << std::endl;
  saveGeomechDataNewKeywords(output_path + data.config.domain_file);

  if (!data.vEfrac.empty())
  {
    std::cout << "save embedded fractures" << std::endl;
    saveEmbeddedFractureProperties(output_path + data.config.efrac_file);
  }

  std::cout << "save mech boundary conditions" << std::endl;
  saveBoundaryConditions(output_path + data.config.bcond_file);

  if (data.dfm_faces.size() > 0)
  {
    std::cout << "save discrete fractures" << std::endl;
    saveDiscreteFractureProperties(output_path + data.config.discrete_frac_file);
  }

  std::cout << "save flow data" << std::endl;
  CalcTranses::save_output(data.flow_data, output_path);
}


void OutputData::saveGeometry(const std::string & output_path)
{
  string outstring;

  ofstream geomechfile;

  // GEOMETRY
  outstring =   output_path + "gm_geometry.txt";
  geomechfile.open(outstring.c_str());
  geomechfile << "GMDIMS" << endl;

  geomechfile << grid.n_vertices() << "\t"
              << grid.n_cells() << "\t"
              << grid.n_faces();
  geomechfile << "/" << endl << endl;

  const std::string vtk_file = output_path + "reservoir_mesh.vtk";

  IO::VTKWriter::write_vtk(grid.vertices.points, grid.cells,
                           grid.shape_ids, vtk_file);

  if (data.dfm_faces.size() > 0)
  { // DFM frac geometry
    mesh::SurfaceMesh<double> frac_msh(/* tol = */ 1e-6);
    for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
      if(data.is_fracture(face.marker()))  // dfm frac
      {
        angem::Polygon<double> poly_face(grid.vertices.points,
                                         face.vertex_indices());
        frac_msh.insert(poly_face);
      }

    const std::string vtk_dfm_file = output_path + "dfm.vtk";
    IO::VTKWriter::write_vtk(frac_msh.vertices.points, frac_msh.polygons,
                             vtk_dfm_file);

  }

  geomechfile.precision(12);
  cout << "write all coordinates\n";
  geomechfile << "GMNODE_COORDS" << endl;
  for (const auto & vertex : grid.vertices)
      geomechfile << vertex[0] << "\t"
                  << vertex[1] << "\t"
                  << vertex[2] << "\n";

  cout << "write all elements\n";
  geomechfile << "GMCELL_NODES" << endl;
  for (auto cell=grid.begin_cells(); cell!=grid.end_cells(); ++cell)
  {
    const auto & vertices = cell.vertices();
    geomechfile << vertices.size() << "\t";

    switch (cell.shape_id())
    {
      case 25: // super wierd element 25
        {
          for (int j = 0; j < 8; j++)
            geomechfile << vertices[j] + 1 << "\t";
          geomechfile << vertices[8] + 1 << "\t";
          geomechfile << vertices[11] + 1 << "\t";
          geomechfile << vertices[13] + 1 << "\t";
          geomechfile << vertices[9] + 1 << "\t";

          geomechfile << vertices[16] + 1 << "\t";
          geomechfile << vertices[18] + 1 << "\t";
          geomechfile << vertices[19] + 1 << "\t";
          geomechfile << vertices[17] + 1 << "\t";

          geomechfile << vertices[10] + 1 << "\t";
          geomechfile << vertices[12] + 1 << "\t";
          geomechfile << vertices[14] + 1 << "\t";
          geomechfile << vertices[15] + 1 << "\t";
          break;
        }
      case 26:
        {
          for (int j = 0; j < 6; j++)
            geomechfile << vertices[j] + 1 << "\t";

          geomechfile << vertices[6] + 1 << "\t";
          geomechfile << vertices[9] + 1 << "\t";
          geomechfile << vertices[7] + 1 << "\t";

          geomechfile << vertices[12] + 1 << "\t";
          geomechfile << vertices[14] + 1 << "\t";
          geomechfile << vertices[13] + 1 << "\t";

          geomechfile << vertices[8] + 1 << "\t";
          geomechfile << vertices[10] + 1 << "\t";
          geomechfile << vertices[11] + 1 << "\t";
          break;
        }
      default:
        {
          for (const auto vertex : vertices)
            geomechfile << vertex + 1 << "\t";
          break;
        }
    }

    geomechfile << std::endl;
  }

  geomechfile << "/" << endl << endl;

  geomechfile << "GMCELL_TYPE" << endl;
  for (auto cell=grid.begin_cells(); cell!=grid.end_cells(); ++cell)
    geomechfile <<  cell.shape_id() << std::endl;
  geomechfile << "/" << std::endl << std::endl;

  geomechfile << "GMCELL_TO_FLOWCELLS" << endl;
  for (const auto & flow_cells : data.gm_cell_to_flow_cell)
  {
    const std::size_t n_connected_elements = flow_cells.size();
    if (n_connected_elements == 0)
      geomechfile << 1 << "\t" << -1 << std::endl;
    else
    {
      geomechfile << n_connected_elements << "\t";
      for (const std::size_t ielement : flow_cells)
        geomechfile << ielement + 1 << "\t";
      geomechfile << std::endl;
    }
  }
  geomechfile << "/\n\n";

  std::cout << "write all faces\n";

  geomechfile << "GMFACE_NODES\n";
  for (auto face=grid.begin_faces(); face!=grid.end_faces(); ++face)
  {
    const std::vector<std::size_t> ivertices = face.vertex_indices();
    geomechfile << ivertices.size() << "\t";
    for (const std::size_t ivertex : ivertices)
      geomechfile << ivertex + 1 << "\t";
    geomechfile << std::endl;
  }

  geomechfile << "/" << std::endl << std::endl;

  geomechfile << "GMFACE_TYPE" << std::endl;
  for (auto face=grid.begin_faces(); face!=grid.end_faces(); ++face)
    if (face.vtk_id() != -1)
      geomechfile << face.vtk_id() << std::endl;
  geomechfile << "/" << std::endl << std::endl;

  geomechfile << "GMFACE_GMCELLS" << std::endl;
  for (auto face=grid.begin_faces(); face!=grid.end_faces(); ++face)
  {
    const auto & neighbors = face.neighbors();
    geomechfile << neighbors.size() << std::endl;
    for (const auto & neighbor : neighbors)
      geomechfile << neighbor + 1 << std::endl;
  }
  geomechfile << "/" << std::endl << std::endl;

    // write vtk data
    // std::cout  << "Writing SDA props" << std::endl;
  if (!data.vEfrac.empty())
  {
    std::size_t n_efrac_vertices = 0;
    // make up a vector of all sda vertices
    for (const auto & efrac : data.vEfrac)
      n_efrac_vertices += efrac.mesh.vertices.size();

    std::vector<angem::Point<3,double>> efrac_verts(n_efrac_vertices);
    std::size_t ivertex = 0;
    for (const auto & efrac : data.vEfrac)
      for (const auto & p : efrac.mesh.vertices)
      {
        efrac_verts[ivertex] = p;
        ivertex++;
      }

    std::size_t n_efrac_elements = 0;
    std::size_t vind_size_total = 0;

    for (const auto & efrac : data.vEfrac)
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
    for (const auto & efrac : data.vEfrac)
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
}


void OutputData::saveGeomechDataNewKeywords(const std::string file_name)
{
  std::cout << "Writing all domain from user input properties" << std::endl;
  {  // write domain properties
    std::ofstream geomechfile;
    geomechfile.open(file_name.c_str());

    for (std::size_t ivar=0; ivar<data.rockPropNames.size(); ++ivar)
    {
      if ( data.config.expression_type[ivar] != 1 )  // only mechanics kwds
        continue;

      geomechfile << data.rockPropNames[ivar] << std::endl;
      for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
      {
        const std::size_t icell = cell.index();
        geomechfile << data.vsCellRockProps[icell].v_props[ivar] << "\t";
        if ((icell + 1) % n_entries_per_line == 0)
          geomechfile << std::endl;
      }
      geomechfile << std::endl << "/" << std::endl << std::endl;
    }

    geomechfile.close();
  }



  // //   geomechfile <<  "COMPDAT\n";
  // //   for ( int iw = 0; iw < pSim->nWells; iw++ )
  // //   {
  // //     for ( int i = 0; i < pSim->vsWell[iw].vID.size(); ++i )
  // //       geomechfile <<  "W" << iw << "\t" << pSim->vsWell[iw].vID[i]+1 <<" 1 1 1 OPEN 1* " << pSim->vsWell[iw].vWi[i] << " 4* Z/" << endl;
  // //   }
  // //   geomechfile << "/\n\n";
  // //   geomechfile.close();
  // // }
}

void OutputData::saveEmbeddedFractureProperties(const std::string file_name)
{
  std::cout  << "Writing SDA props" << std::endl;
  std::ofstream geomechfile;
  geomechfile.open(file_name.c_str());

  geomechfile << "GM_EFRAC_CELLS" << std::endl;
  for (const auto & efrac : data.vEfrac)
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
  for (const auto & efrac : data.vEfrac)
    for (std::size_t i=0; i<efrac.points.size(); ++i)
      geomechfile << efrac.points[i] << std::endl;
  geomechfile << "/" << std::endl << std::endl;

  geomechfile << "GM_EFRAC_DIP" << std::endl;
  for (const auto & efrac : data.vEfrac)
    for (std::size_t i=0; i<efrac.points.size(); ++i)
    {
      geomechfile << efrac.dip[i] << "\t";
      if ((i+1) % n_entries_per_line == 0)
        geomechfile << std::endl;
    }
  geomechfile << "/" << std::endl << std::endl;

  geomechfile << "GM_EFRAC_STRIKE" << std::endl;
  for (const auto & efrac : data.vEfrac)
    for (std::size_t i=0; i<efrac.points.size(); ++i)
    {
      geomechfile << efrac.strike[i] << "\t";
      if ((i+1)%n_entries_per_line == 0) geomechfile << std::endl;
    }
  geomechfile << "/" << std::endl << std::endl;

  geomechfile << "GM_EFRAC_COHESION" << std::endl;
  for (const auto & efrac : data.vEfrac)
    geomechfile << efrac.cohesion << std::endl;
  geomechfile << "/" << std::endl << std::endl;

  geomechfile << "GM_EFRAC_FRICTION" << std::endl;
  for (const auto & efrac : data.vEfrac)
    geomechfile << efrac.friction_angle << std::endl;
  geomechfile << "/" << std::endl << std::endl;

  geomechfile << "GM_EFRAC_DILATION" << std::endl;
  for (const auto & efrac : data.vEfrac)
    geomechfile << efrac.dilation_angle << std::endl;
  geomechfile << "/" << std::endl << std::endl;

  geomechfile.close();
}


void OutputData::saveBoundaryConditions(const std::string file_name)
{
  std::ofstream geomechfile;
  std::cout << "Computing Dirichlet nodes" << std::endl;
  // /* this chunk of code find dirichlet nodes
  //  * Algorithm:
  //  * 1. First add all Dirichlet nodes specified by the user
  //  * to the set
  //  * 2. loop Dirichlet faces.
  //  * loop nodes withing dirichlet faces
  //  * add node as a part of dirichlet face
  //  */
  const int dim = 3;
  std::vector<std::unordered_set<std::size_t>> setDisp(dim);
  std::vector<std::vector<std::size_t>> vv_disp_node_ind(dim);
  std::vector<std::vector<double>> vv_disp_node_values(dim);

  // search tolerance
  const double tol = data.config.node_search_tolerance;

  for (std::size_t ivert=0; ivert<grid.vertices.size(); ++ivert)
  {
    const auto & vertex = grid.vertices[ivert];

    for (const auto & bc_node : data.config.bc_nodes)
    {
      if (bc_node.coord.distance(vertex) < tol)
      {
        for (int d=0; d<dim; ++d)
          if (bc_node.value[d] != data.config.nan)
          {
            setDisp[d].insert(ivert);
            vv_disp_node_values[d].push_back(bc_node.value[d]);
            vv_disp_node_ind[d].push_back(ivert);
          }
      }
    }
  }

  if ( data.nDirichletFaces > 0 )
    for (auto face=grid.begin_faces(); face!=grid.end_faces(); ++face)
      if (data.is_boundary(face.marker()))
      {
        auto facet_it = data.boundary_faces.find(face.index());
        if (facet_it != data.boundary_faces.end())
        {
          if (facet_it->second.ntype == 1) // dirichlet
            for (std::size_t d=0; d<dim; ++d)
              if (facet_it->second.condition[d] != data.config.nan)
                for (const auto ivertex : face.vertex_indices())
                {
                  const auto ret = setDisp[d].insert(ivertex);
                  // check if already in set
                  if(ret.second)
                  {
                    vv_disp_node_ind[d].push_back(ivertex);
                    vv_disp_node_values[d].push_back(facet_it->second.condition[d]);
                  } // end if vertex is new
                }
        }
        else  // should not happen
        {
          std::cout << "face " << face.index()
                    << " is not a boundary face!! Aborting!!!"
                    << std::endl;
          abort();
        }
      }

  setDisp.clear();

  geomechfile.open(file_name.c_str());

  // save dirichlet faces
  for (std::size_t j=0; j<dim; ++j)
    if (vv_disp_node_ind[j].size() > 0)
    {
      std::cout << "write all Dirichlet nodes" << std::endl;
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

  if ( data.nNeumannFaces > 0 )
  {
    cout << "write all Neumann faces\n";

    geomechfile << "GMFACE_TRACTION_TXYZ" << endl;
    // for ( std::size_t i = 0; i < data.nPhysicalFacets; i++ )
    for (const auto facet_it : data.boundary_faces)
    {
      const auto & facet = facet_it.second;
      if (facet.ntype == 2)  // neumann
      {
        geomechfile << facet.nface + 1 << "\t";
        geomechfile << facet.condition << std::endl;
      }
    }
    geomechfile << "/\n\n";
  }

  geomechfile.close();
}


void OutputData::saveDiscreteFractureProperties(const std::string file_name)
{
  std::cout << "write discrete fracs" << std::endl;

  std::ofstream geomechfile;
  geomechfile.open(file_name.c_str());
  set<int>::iterator itsetint;

  int counter = 0;
  int nFractures_ = 0;
  cout << "write all fractured faces\n";

  geomechfile << "GMFACE_FRACTURE_TO_FLOWCELL\n";
  for (const auto facet_it : data.dfm_faces)
  {
    geomechfile << facet_it.second.nface + 1 << "\t";
    geomechfile << facet_it.second.nfluid + 1 << endl;
  }
  // for(itsetint = pSim->setIdenticalInternalMarker.begin();
  //     itsetint != pSim->setIdenticalInternalMarker.end();
  //     itsetint++, nFractures_++)
  // {
  //   for (int i = 0; i < pSim->nPhysicalFacets; i++)
  //   {
  //     if( pSim->vsPhysicalFacet[i].nmark == *itsetint )
  //     {
  //       geomechfile << pSim->vsPhysicalFacet[i].nface + 1 << "\t";
  //       geomechfile << pSim->vsPhysicalFacet[i].nfluid + 1 << endl;
  //       if( pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].nNeighbors !=2 )
  //       {
  //         cout << "Fracture interface # " << nFractures_ << endl;
  //         cout << "Global interface   # " << pSim->vsPhysicalFacet[i].nface << endl;
  //         cout << "Number od neighbors  " << pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].nNeighbors << endl;
  //         cout << "Wrong msh file. Mesh verticies are not connected on fracture interface" << endl;
  //         exit(0);
  //       }
  //     }
  //   }
  // }
  geomechfile << "/" << endl << endl;

  //   nFractures_ = 0;
  geomechfile << "GMFACE_FRACTURE_CONDUCTIVITY" << std::endl;
  for (const auto facet_it : data.dfm_faces)
    geomechfile << facet_it.second.conductivity << std::endl;

  geomechfile << "GMFACE_FRACTURE_REGION" << std::endl;
  for (const auto facet_it : data.dfm_faces)
    geomechfile << 1 << std::endl;
  geomechfile << "/" << std::endl << std::endl;

  geomechfile << "GMFACE_FRACTURE_GROUP" << std::endl;
  for (const auto facet_it : data.dfm_faces)
    geomechfile << 1 << std::endl;
  geomechfile << "/" << std::endl << std::endl;

  geomechfile.close();
}
