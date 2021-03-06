#include "OutputDataGPRS.hpp"
#include "VTKWriter.hpp"

#include <sys/stat.h>

namespace gprs_data
{

const int n_entries_per_line = 10;


OutputDataGPRS::OutputDataGPRS(SimData & sim_data,
                       mesh::Mesh & grid)
    :
    data(sim_data),
    grid(grid),
    ordered_faces(grid.get_ordered_faces())
{}


void OutputDataGPRS::write_output(const std::string & output_path)
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

  std::cout << "save mech boundary conditions: "
            << output_path + data.config.bcond_file
            << std::endl;
  saveBoundaryConditions(output_path + data.config.bcond_file);

  if (data.dfm_faces.size() > 0)
  {
    std::cout << "save discrete fractures" << std::endl;
    saveDiscreteFractureProperties(output_path + data.config.discrete_frac_file);
  }

  if (!data.wells.empty())
  {
    std::cout << "save wells" << std::endl;
    saveWells(output_path + data.config.wells_file);
  }

  // flow discretization
  std::cout << "save flow discretization" << std::endl;
  flow::CalcTranses::save_output(data.flow_data, output_path);

  // multiscale
  if (data.ms_flow_data.partitioning.size() > 0)
    saveFlowMultiScaleData(output_path + data.config.flow_ms_file);
  if (data.ms_mech_data.partitioning.size() > 0)
    saveMechMultiScaleData(output_path + data.config.mech_ms_file);
}


void OutputDataGPRS::saveGeometry(const std::string & output_path)
{
  // gprs output

  // GEOMETRY
  const std::string outstring = output_path + data.config.mechanics_domain_file;
  std::cout << "writing file " << outstring << std::endl;

  std::ofstream geomechfile;
  geomechfile.open(outstring.c_str());
  geomechfile << "GMDIMS" << endl;

  geomechfile << grid.n_vertices() << "\t"
              << grid.n_cells() << "\t"
              << grid.n_faces() - data.dfm_faces.size();  // not counting the split faces
  geomechfile << "/" << endl << endl;

  geomechfile.precision(6);
  cout << "write all coordinates\n";
  geomechfile << "GMNODE_COORDS" << endl;
  for (const auto & vertex : grid.vertices)
      geomechfile << vertex[0] << "\t"
                  << vertex[1] << "\t"
                  << vertex[2] << "\n";
  geomechfile << "/" << std::endl << std::endl;

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
  for (const auto & face : ordered_faces)
  {
    if (data.is_fracture(face.marker()))  // skip non-master faces
      if (face.index() != face.master_index())
        continue;

    const std::vector<std::size_t> ivertices = face.vertex_indices();
    geomechfile << ivertices.size() << "\t";
    for (const std::size_t ivertex : ivertices)
      geomechfile << ivertex + 1 << "\t";
    geomechfile << std::endl;
  }

  geomechfile << "/" << std::endl << std::endl;

  geomechfile << "GMFACE_TYPE" << std::endl;
  for (const auto & face : ordered_faces)
  {
    if (data.is_fracture(face.marker()))  // skip non-master faces
      if (face.index() != face.master_index())
        continue;
    geomechfile << face.vtk_id() << std::endl;
  }
  geomechfile << "/" << std::endl << std::endl;

  std::cout << "writing face-cell connection" << std::endl;
  geomechfile << "GMFACE_GMCELLS" << std::endl;
  for (const auto & face : ordered_faces)
  {
    if (data.is_fracture(face.marker()))  // timur want to retain neighbors of master frac face
    {
      if (face.index() == face.master_index())
      {
        const auto it_frac_face = data.dfm_faces.find(face.master_index());
        if (it_frac_face == data.dfm_faces.end())
        {
          std::cout << "bug in dfm connections" << std::endl;
          exit(0);
        }
        const auto & neighbors = it_frac_face->second.neighbor_cells;

        geomechfile << neighbors.size() << "\t";
        for (const auto & neighbor : neighbors)
          geomechfile << neighbor + 1 << "\t";
        geomechfile << std::endl;
      }
    }
    else
    {
      assert( face.index() == face.master_index() );
      const auto & neighbors = face.neighbors();
      geomechfile << neighbors.size() << "\t";
      for (const std::size_t neighbor : neighbors)
        geomechfile << neighbor + 1 << "\t";
      geomechfile << std::endl;
    }
  }
  geomechfile << "/" << std::endl << std::endl;
  geomechfile.close();
}


void OutputDataGPRS::saveGeomechDataNewKeywords(const std::string file_name)
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
}


void OutputDataGPRS::saveEmbeddedFractureProperties(const std::string file_name)
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


void OutputDataGPRS::saveBoundaryConditions(const std::string file_name)
{
  std::ofstream geomechfile;
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

  if ( data.n_dirichlet_faces > 0 )
    for (const auto & face : ordered_faces)
      if (data.is_boundary(face.marker()))
      {
        const auto facet_it = data.boundary_faces.find(face.master_index());
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

      for(std::size_t i = 0; i < vv_disp_node_ind[j].size(); ++i)
      {
        geomechfile << vv_disp_node_ind[j][i] + 1 << "\t";
        geomechfile << vv_disp_node_values[j][i] << std::endl;
      }
      geomechfile << "/" << std::endl << std::endl;
    }

  if ( data.n_neumann_faces > 0 )
  {
    std::cout << "write all Neumann faces" << std::endl;

    geomechfile << "GMFACE_TRACTION_TXYZ" << std::endl;
    for (const auto & facet_it : data.boundary_faces)
      if (facet_it.second.ntype == 2)  // neumann
      {
        geomechfile << facet_it.second.nface + 1 << "\t";
        geomechfile << facet_it.second.condition << std::endl;
      }

    geomechfile << "/" << std::endl << std::endl;
  }

  geomechfile.close();
}


void OutputDataGPRS::saveDiscreteFractureProperties(const std::string file_name)
{
  std::cout << "write discrete fracs" << std::endl;

  std::ofstream geomechfile;
  geomechfile.open(file_name.c_str());
  set<int>::iterator itsetint;

  cout << "write all fractured faces\n";
  geomechfile << "GMFACE_FRACTURE_TO_FLOWCELL" << std::endl;
  for (const auto & face_it : data.dfm_faces)
  {
    geomechfile << face_it.second.nface + 1 << "\t";
    if (face_it.second.coupled)
      geomechfile << face_it.second.nfluid + 1 << std::endl;
    else
      geomechfile << -1 << std::endl;
  }
  geomechfile << "/" << std::endl << std::endl;

  geomechfile << "GMFACE_FRACTURE_CONDUCTIVITY" << std::endl;
  for (const auto facet_it : data.dfm_faces)
    geomechfile << facet_it.second.conductivity << std::endl;
  geomechfile << "/" << std::endl << std::endl;

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


void OutputDataGPRS::saveWells(const std::string file_name)
{
  std::ofstream file;
  file.open(file_name.c_str());

  file << "WELSPECS" << std::endl;
  for (const auto & well : data.wells)
  {
    if (well.connected_volumes.empty())
    {
      std::cout << "WARNING. The well " << well.name
                << " is not connected to any CVs" << std::endl;
      continue;
    }

    file << well.name << "\tGROUP1\t";
    // connected volume + j + k empty
    file << well.connected_volumes[0] << " 1 1 ";
    // reference depth
    file << -well.reference_depth << " /" << std::endl;
  }
  file << "/" << std::endl << std::endl;

  file << "COMPDAT" << std::endl;
  for (const auto & well : data.wells)
  {
    for (std::size_t i=0; i<well.connected_volumes.size(); ++i)
    {
      file << well.name << "\t";
      file << well.connected_volumes[i] + 1 << "\t";
      // j, k1:k2 open sat_table_number
      file << "1\t1\t1\tOPEN\t1*\t";
      file << well.indices[i] * flow::CalcTranses::transmissibility_conversion_factor << "\t";
      file << 2*well.radius << "\t";
      file << "/" << std::endl;
    }
  }
  file << "/" << std::endl << std::endl;

  file.close();
}


void OutputDataGPRS::saveFlowMultiScaleData(const std::string file_name)
{

}


void OutputDataGPRS::saveMechMultiScaleData(const std::string file_name)
{
  std::ofstream out;
  out.open(file_name.c_str());
  const auto & ms = data.ms_mech_data;

  // save partitioing
  out << "GMMSPARTITIONING";
  for (std::size_t i=0; i < ms.partitioning.size(); ++i)
  {
    if (i % n_entries_per_line == 0) out << endl;
    out << ms.partitioning[i] << " ";
  }
  out << "/" << endl << endl;

  // save support
  out << "GMMSSUPPORT ";
  for (std::size_t i=0; i < ms.n_coarse; ++i)
  {
    out << endl;
    out << ms.support_internal[i].size() << " "  // number of cells (centroid)
        << ms.support_boundary[i].size() << " "; // number of boundary nodes

    // first write the centroid (vertex)
    out << ms.centroids[i] << " ";

    // internal cells
    size_t counter = 3;
    for (const size_t cell : ms.support_internal[i])
    {
      if (counter++ % n_entries_per_line == 0) out << endl;
      out << cell << " ";
    }

    // boundary nodes
    for (const size_t vertex : ms.support_boundary[i])
    {
      if (counter++ % n_entries_per_line == 0) out << endl;
      out << vertex << " ";
    }
  }

  out << "/" << endl << endl;

  out.close();
}

}
