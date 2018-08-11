#include "femout.hpp"

OutputData::OutputData(SimData * pSimData)
{
  pSim = pSimData;
}

OutputData::~OutputData()
{
}



#include <sys/stat.h>
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
  outstring =   output_path + "fl_dimens.txt";
  geomechfile.open(outstring.c_str());

  geomechfile << "DIMENS" << endl;
  geomechfile << pSim->nInternalBoundaryFaces + pSim->nCells << "\t" << 1 << "\t" << 1 << " /" << endl;
  geomechfile.close();


  outstring = output_path + "gm_depth.txt";
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

    const std::size_t shift = pSim->n_default_vars();
    for (std::size_t ivar=0; ivar<pSim->config.all_vars.size() - shift; ++ivar)
    {
      geomechfile << pSim->config.all_vars[shift+ivar] << std::endl;
      for (std::size_t icell=0; icell<pSim->nCells; ++icell)
      {
        geomechfile << pSim->vsCellRockProps[icell].v_props[ivar] << "\t";
        if ((icell+1)%10 == 0)
          std::cout << std::endl;
      }
      geomechfile << std::endl << "/" << std::endl << std::endl;
    }

    geomechfile.close();
  }

  {  // Strong discontinuity
    std::cout  << "Writing SDA props" << std::endl;
    outstring =   output_path + "gm_SDA.txt";
    geomechfile.open(outstring.c_str());

    std::size_t ef_ind = 0;
    const auto & efrac = pSim->vsEmbeddedFractures[ef_ind];
    const std::size_t n_sda = efrac.cells.size();

    { // cells
      geomechfile << "GM_EFRAC_CELLS\n";
      geomechfile << n_sda << std::endl;

      for (std::size_t i=0; i<n_sda; ++i)
      {
        geomechfile << efrac.cells[i] + 1 << "\t";
        if ((i+1)%10 == 0)
          geomechfile << std::endl;
      }
      geomechfile << "/" << std::endl << std::endl;
    }
    { // points
      const std::size_t dim = 3;
      geomechfile << "GM_EFRAC_POINTS" << std::endl;
      for (std::size_t i=0; i<n_sda; ++i)
      {
        for (std::size_t j=0; j<dim; ++j)
          geomechfile << efrac.points[i](j) << "\t";
        geomechfile << std::endl;
      }
      geomechfile << "/" << std::endl << std::endl;
    }
    { // Dip
      geomechfile << "GM_EFRAC_DIP" << std::endl;
      for (std::size_t i=0; i<n_sda; ++i)
      {
        geomechfile << efrac.dip[i] << "\t";
        if ((i+1)%10 == 0)
          geomechfile << std::endl;
      }
      geomechfile << "/" << std::endl << std::endl;
    }
    { // Strike
      geomechfile << "GM_EFRAC_STRIKE" << std::endl;
      for (std::size_t i=0; i<n_sda; ++i)
      {
        geomechfile << efrac.strike[i] << "\t";
        if ((i+1)%10 == 0)
          geomechfile << std::endl;
      }
      geomechfile << "/" << std::endl << std::endl;
    }
    {// cohesion
      geomechfile << "GM_EFRAC_COHESION" << std::endl;
      geomechfile << efrac.cohesion << std::endl;
      geomechfile << "/" << std::endl << std::endl;
    }
    {// friction
      geomechfile << "GM_EFRAC_FRICTION" << std::endl;
      geomechfile << efrac.friction_angle << std::endl;
      geomechfile << "/" << std::endl << std::endl;
    }
    {// dilation
      geomechfile << "GM_EFRAC_DILATION" << std::endl;
      geomechfile << efrac.dilation_angle << std::endl;
      geomechfile << "/" << std::endl << std::endl;
    }

    geomechfile.close();
  }

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


  // Dirichlet faces
  // @HACK SUPPORT POINTS
  std::set<int>::iterator it;
  const int dim = 3;
  std::vector<set<int>> setDisp(dim);
  std::vector<std::vector<std::size_t>> vv_disp_node_ind(dim);
  std::vector<std::vector<double>> vv_disp_node_values(dim);

  setDisp.clear();
  if ( pSim->nDirichletFaces > 0 )
    for ( int i = 0; i < pSim->nPhysicalFacets; i++ )
      if ( pSim->vsPhysicalFacet[i].ntype == 1 )  // dirichlet
        for (std::size_t j=0; j<dim; ++j)
        {
          if ( pSim->vsPhysicalFacet[i].condition[j] != pSim->dNotNumber )
          {
            int nvrtx = pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].nVertices;
            for ( int ivrtx = 0; ivrtx < nvrtx; ++ivrtx )
            {
              // check if already in set
              const auto ret = setDisp[j].insert
                  (pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].vVertices[ivrtx]);
              if(ret.second)
              {
                std::size_t vertex_ = pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].vVertices[ivrtx];
                vv_disp_node_ind[j].push_back(vertex_);
                vv_disp_node_values[j].push_back(pSim->vsPhysicalFacet[i].condition[j]);
              }
            }  // end face vertices loop
          } // end nan condition
        }  // end j loop

  setDisp.clear();

  cout << "write all Dirichlet faces\n";
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

  // if(vDisp_x_Idx.size() > 0)
  // {
  //   std::cout << "gmnode_bdispx" << std::endl << std::flush;

  //   geomechfile << "GMNODE_BCDISPX\n";
  //   for(int i = 0; i < vDisp_x_Idx.size(); ++i)
  //   {
  //     geomechfile <<  vDisp_x_Idx[i] + 1 << "\t";
  //     geomechfile << vDisp_x_Val[i] << endl;
  //   }
  //   geomechfile << "/\n\n";
  // }
  // if(vDisp_y_Idx.size() > 0)
  // {
  //   std::cout << "gmnode_bdispy" << std::endl << std::flush;
  //   geomechfile << "GMNODE_BCDISPY\n";
  //   for(int i = 0; i < vDisp_y_Idx.size(); ++i)
  //   {
  //     geomechfile <<  vDisp_y_Idx[i] + 1 << "\t";
  //     geomechfile << vDisp_y_Val[i] << endl;
  //   }
  //   geomechfile << "/\n\n";
  // }

  // if(vDisp_z_Idx.size() > 0)
  // {
  //   geomechfile << "GMNODE_BCDISPZ\n";
  //   for(int i = 0; i < vDisp_z_Idx.size(); ++i)
  //   {
  //     geomechfile <<  vDisp_z_Idx[i] + 1 << "\t";
  //     geomechfile << vDisp_z_Val[i] << endl;
  //   }
  //   geomechfile << "/\n\n";
  // }
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

  std::cout << "pressure" << std::endl;

  outstring =   output_path + "fl_pres.txt";
  geomechfile.open(outstring.c_str());

  // fractures first
  geomechfile << "PRESSURE" << endl;
  for(int iface = 0; iface < pSim->nFaces; iface++)
  {
    if( pSim->vsFaceCustom[iface].nMarker > 0)
      geomechfile << (pSim->vsCellRockProps[pSim->vsFaceCustom[iface].vNeighbors[0]].pressure + pSim->vsCellRockProps[pSim->vsFaceCustom[iface].vNeighbors[1]].pressure) / 2.0 << endl;
  }
  // matrix
  for(int ib = 0; ib < pSim->nCells; ++ib)
    geomechfile << pSim->vsCellRockProps[ib].pressure << endl;
  geomechfile << "\n/\n\n";
  geomechfile.close();

  std::cout << "temperature" << std::endl;

  outstring =   output_path + "fl_temp.txt";
  geomechfile.open(outstring.c_str());

  // fractures first
  geomechfile << "RTEMP" << endl;
  for(int iface = 0; iface < pSim->nFaces; iface++)
  {
    // internal suface
    if( pSim->vsFaceCustom[iface].nMarker > 0)
    {
      int n1_ = pSim->vsFaceCustom[iface].vNeighbors[0];
      int n2_ = pSim->vsFaceCustom[iface].vNeighbors[1];
      geomechfile << min(pSim->vsCellRockProps[n1_].temp, pSim->vsCellRockProps[n2_].temp) << endl;
    }
  }
  // matrix
  for(int ib = 0; ib < pSim->nCells; ++ib)
    geomechfile << pSim->vsCellRockProps[ib].temp << endl;

  geomechfile << "\n/\n\n";
  geomechfile.close();

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


  std::cout << "saturations" << std::endl;

  outstring =   output_path + "fl_satnum.txt";
    geomechfile.open(outstring.c_str());

  geomechfile << "SATNUM" << endl;
  for(int iface = 0; iface < pSim->nFaces; iface++)
  {
    if( pSim->vsFaceCustom[iface].nMarker > 0)
      geomechfile << 1 << endl;
  }

  for(int ib = 0; ib < pSim->nCells; ++ib)
    geomechfile << 0 << endl;

  geomechfile << "/" << endl;
  geomechfile.close();


  // outstring =   output_path + "fl_thc.txt";
  //   geomechfile.open(outstring.c_str());

  // geomechfile << "THCROCK" << endl;
  // for(int iface = 0; iface < pSim->nFaces; iface++)
  // {
  //   // internal suface
  //   if( pSim->vsFaceCustom[iface].nMarker > 0)
  //   {
  //     int n1_ = pSim->vsFaceCustom[iface].vNeighbors[0];
  //     int n2_ = pSim->vsFaceCustom[iface].vNeighbors[1];
  //     geomechfile << min(pSim->vsCellRockProps[n1_].thc, pSim->vsCellRockProps[n2_].thc) << endl;
  //   }
  // }
  // for(int ib = 0; ib < pSim->nCells; ++ib)
  //   geomechfile << pSim->vsCellRockProps[ib].thc << endl;

  // geomechfile << "/" << endl;
  // geomechfile.close();


  cout << "write all wells\n";
  if ( pSim->nWells > 0 )
  {
    outstring =   output_path + "fl_wells.txt";
    geomechfile.open ( outstring.c_str() );
    geomechfile <<  "WELSPECS\n";

    for ( int iw = 0; iw < pSim->nWells; iw++ )
    {
      if(pSim->vsWell[iw].vID.size() > 0)
	geomechfile <<  "W" << iw << " 1* " << pSim->vsWell[iw].vID[0]+1 << " 1 " <<  pSim->vsWell[iw].datum << " WATER /" << endl;
    }

    geomechfile << "/\n\n";

    geomechfile <<  "COMPDAT\n";
    for ( int iw = 0; iw < pSim->nWells; iw++ )
    {
      for ( int i = 0; i < pSim->vsWell[iw].vID.size(); ++i )
        geomechfile <<  "W" << iw << "\t" << pSim->vsWell[iw].vID[i]+1 <<" 1 1 1 OPEN 1* " << pSim->vsWell[iw].vWi[i] << " 4* Z/" << endl;
    }
    geomechfile << "/\n\n";
    geomechfile.close();
  }
}
