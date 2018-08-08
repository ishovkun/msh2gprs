#include "femout.hpp"

OutputData::OutputData(SimData * pSimData)
{
  pSim = pSimData;
}

OutputData::~OutputData()
{
}



#include <sys/stat.h>
void OutputData::writeGeomechDataNewKeywords()
{
  stringstream out;
  out << "model";
  out << "/";
  string outputPath_ = ""; // out.str();

#if 0 //defined(_WIN32)
  _mkdir(outputPath_.c_str()); // create directory
//#else
  mkdir(outputPath_.c_str(), 0777);
#endif
  string outstring = pSim->outstream;

  ofstream geomechfile;
  outstring =   outputPath_ + "fl_dimens.txt";
  geomechfile.open(outstring.c_str());

  geomechfile << "DIMENS" << endl;
  geomechfile << pSim->nInternalBoundaryFaces + pSim->nCells << "\t" << 1 << "\t" << 1 << " /" << endl;
  geomechfile.close();


  outstring = outputPath_ + "gm_depth.txt";
  geomechfile.open(outstring.c_str());
  geomechfile << "GDEPTH" << endl;

  for (int i = 0; i < pSim->nCells; i++)
  {
	  geomechfile << pSim->vsCellCustom[i].vCenter [2] << "\n";
  }
  geomechfile << "/" << endl << endl;
  geomechfile.close();

  // GEOMETRY
  outstring =   outputPath_ + "gm_geometry.txt";
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

 outstring =   outputPath_ + "gm_model.txt";
 geomechfile.open(outstring.c_str());
  geomechfile << "GMCELL_MODEL" << endl;
  for(int ib = 0; ib < pSim->nCells; ++ib)
  {
    geomechfile << pSim->vsCellRockProps[ib].model << endl;
  }
  geomechfile << "/" << endl << endl;
 geomechfile.close();

 outstring =   outputPath_ + "gm_capacity.txt";
 geomechfile.open(outstring.c_str());
  geomechfile << "GMCELL_HEAT_CAPACITY" << endl;
  for(int ib = 0; ib < pSim->nCells; ++ib)
  {
    geomechfile << pSim->vsCellRockProps[ib].heat_capacity << endl;
  }
  geomechfile << "/" << endl << endl;
 geomechfile.close();

 outstring =   outputPath_ + "gm_density.txt";
 geomechfile.open(outstring.c_str());

  geomechfile << "GMCELL_DENSITY\n";
  int j = 0;
  for ( int i = 0; i < pSim->nCells; i++, j++ )
  {
    geomechfile << pSim->vsCellRockProps[i].density << "\t";
    if( j == 10 )
    {
      j = 0;
      geomechfile << endl;
    }
  }
  geomechfile << "\n/\n\n";
  geomechfile.close();

 outstring =   outputPath_ + "gm_biot.txt";
 geomechfile.open(outstring.c_str());
  geomechfile << "GMCELL_BIOT\n";
  j = 0;
  for ( int i = 0; i < pSim->nCells; i++, j++ )
  {
    geomechfile << pSim->vsCellRockProps[i].biot << "\t";
    if( j == 10 )
    {
      j = 0;
      geomechfile << endl;
    }
  }
  geomechfile << "\n/\n\n";

  geomechfile << "GMCELL_BIOT_FLOW\n";
  j = 0;
  for ( int i = 0; i < pSim->nCells; i++, j++ )
  {
    geomechfile << pSim->vsCellRockProps[i].biot_flow << "\t";
    if( j == 10 )
    {
      j = 0;
      geomechfile << endl;
    }
  }
  geomechfile << "\n/\n\n";
  geomechfile.close();

  outstring =   outputPath_ + "gm_elastic.txt";
  geomechfile.open(outstring.c_str());
  geomechfile << "GMCELL_YOUNG\n";
  j = 0;
  for ( int i = 0; i < pSim->nCells; i++, j++ )
  {
    geomechfile << pSim->vsCellRockProps[i].young << "\t";
    if( j == 10 )
    {
      j = 0;
      geomechfile << endl;
    }
  }
  geomechfile << "\n/\n\n";

 geomechfile << "GMCELL_POISSON\n";
  j = 0;
  for ( int i = 0; i < pSim->nCells; i++, j++ )
  {
    geomechfile << pSim->vsCellRockProps[i].poisson << "\t";
    if( j == 8)
    {
      j = 0;
      geomechfile << endl;
    }
  }
  geomechfile << "\n/\n\n";
  geomechfile.close();

 outstring =   outputPath_ + "gm_thermal.txt";
 geomechfile.open(outstring.c_str());

 geomechfile << "GMCELL_THERMAL_EXPANSION\n";
  j = 0;
  for ( int i = 0; i < pSim->nCells; i++, j++ )
  {
    geomechfile << pSim->vsCellRockProps[i].thermal_expansion  << "\t";
    if( j == 10 )
    {
      j = 0;
      geomechfile << endl;
    }
  }
  geomechfile << "\n/\n\n";

 geomechfile << "GMCELL_PORE_THERMAL_EXPANSION\n";
  j = 0;
  for ( int i = 0; i < pSim->nCells; i++, j++ )
  {
   geomechfile << 3.0 * pSim->vsCellRockProps[i].thermal_expansion * (pSim->vsCellRockProps[i].biot - pSim->vsCellRockProps[i].poro) << "\t";
   if( j == 10 )
    {
      j = 0;
      geomechfile << endl;
    }
  }
  geomechfile << "\n/\n\n";
  geomechfile.close();

  {  // Strong discontinuity
    std::cout  << "Writing SDA props" << std::endl;
    outstring =   outputPath_ + "gm_SDA.txt";
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

  outstring =   outputPath_ + "gm_reference.txt";
  geomechfile.open(outstring.c_str());

  geomechfile << "GMREF_TEMPERATURE" << endl;
  for(int ib = 0; ib < pSim->nCells; ++ib)
  {
    geomechfile << pSim->vsCellRockProps[ib].ref_temp << endl;
  }
  geomechfile << "/" << endl << endl;

  geomechfile << "GMREF_PRESSURE" << endl;
  for(int ib = 0; ib < pSim->nCells; ++ib)
  {
    geomechfile << pSim->vsCellRockProps[ib].ref_pres << endl;
  }
  geomechfile << "/" << endl << endl;

  geomechfile << "/" << endl << endl;
  geomechfile.close();

  outstring =   outputPath_ + "gm_bcond.txt";
  geomechfile.open(outstring.c_str());

  if ( pSim->nNeumannFaces > 0 )
  {
    cout << "write all Neumann faces\n";

    bool is_ = false;
    for ( int i = 0; i < pSim->nPhysicalFacets; i++ )
    {
      if ( pSim->vsPhysicalFacet[i].ntype == 2 )
        is_= true;
    }
    if ( is_ )
    {
      geomechfile << "GMFACE_TRACTION_N" << endl;
      for ( int i = 0; i < pSim->nPhysicalFacets; i++ )
      {
        if ( pSim->vsPhysicalFacet[i].ntype == 2 )
        {
          geomechfile << pSim->vsPhysicalFacet[i].nface + 1 << "\t";
          geomechfile << pSim->vsPhysicalFacet[i].vCondition[0] << endl;
        }
      }
    }

    is_ = false;
    for ( int i = 0; i < pSim->nPhysicalFacets; i++ )
    {
      if ( pSim->vsPhysicalFacet[i].ntype == 22 )
        is_= true;
    }
    if ( is_ )
    {
      geomechfile << "GMFACE_TRACTION_TXYZ" << endl;
      for ( int i = 0; i < pSim->nPhysicalFacets; i++ )
      {
        if ( pSim->vsPhysicalFacet[i].ntype == 22)
        {
          geomechfile << pSim->vsPhysicalFacet[i].nface + 1 << "\t";
          geomechfile << pSim->vsPhysicalFacet[i].vCondition[0] << "\t";
          geomechfile << pSim->vsPhysicalFacet[i].vCondition[1] << "\t";
          geomechfile << pSim->vsPhysicalFacet[i].vCondition[2] << endl;
        }
      }
    }
    geomechfile << "/\n\n";
  }



  // @HACK SUPPORT POINTS
  std::set<int>::iterator it;
  std::pair<std::set<int>::iterator,bool> ret;
  set<int> setDisp;
  vector<int> vDisp_x_Idx;
  vector<double> vDisp_x_Val;
  vector<int> vDisp_y_Idx;
  vector<double> vDisp_y_Val;
  vector<int> vDisp_z_Idx;
  vector<double> vDisp_z_Val;

  setDisp.clear();
  if ( pSim->nDirichletFaces > 0 )
  {
    setDisp.clear();
    for ( int i = 0; i < pSim->nPhysicalFacets; i++ )
    {
      if ( pSim->vsPhysicalFacet[i].ntype == 1 )
      {
        int n = pSim->vsPhysicalFacet[i].vCondition.size();

        int j = 0;
        if ( pSim->vsPhysicalFacet[i].vCondition[j] != pSim->dNotNumber )
        {
          int nvrtx = pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].nVertices;
          for ( int ivrtx = 0; ivrtx < nvrtx; ++ivrtx )
          {
	    ret = setDisp.insert(pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].vVertices[ivrtx]);
	    if(ret.second == true)
	    {
	      int vertex_ = pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].vVertices[ivrtx];
	      vDisp_x_Idx.push_back(vertex_);
	      vDisp_x_Val.push_back(0.0);
	    }
          }
        }
      }
    }

    setDisp.clear();
    for ( int i = 0; i < pSim->nPhysicalFacets; i++ )
    {
      if ( pSim->vsPhysicalFacet[i].ntype == 1 )
      {
        int n = pSim->vsPhysicalFacet[i].vCondition.size();

        int j = 1;
        if ( pSim->vsPhysicalFacet[i].vCondition[j] != pSim->dNotNumber )
        {
          int nvrtx = pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].nVertices;
          for ( int ivrtx = 0; ivrtx < nvrtx; ++ivrtx )
          {
	    ret = setDisp.insert(pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].vVertices[ivrtx]);
	    if(ret.second == true)
	    {
	      int vertex_ = pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].vVertices[ivrtx];
	      vDisp_y_Idx.push_back(vertex_);
	      vDisp_y_Val.push_back(0.0);
	    }
          }
        }
      }
    }

    setDisp.clear();
    for ( int i = 0; i < pSim->nPhysicalFacets; i++ )
    {
      if ( pSim->vsPhysicalFacet[i].ntype == 1 )
      {
        int n = pSim->vsPhysicalFacet[i].vCondition.size();

        int j = 2;
        if ( pSim->vsPhysicalFacet[i].vCondition[j] != pSim->dNotNumber )
        {
          int nvrtx = pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].nVertices;
          for ( int ivrtx = 0; ivrtx < nvrtx; ++ivrtx )
          {
	    ret = setDisp.insert(pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].vVertices[ivrtx]);
	    int vertex_ = pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].vVertices[ivrtx];
	    if(ret.second == true)
	    {
	      vDisp_z_Idx.push_back(vertex_);
	      vDisp_z_Val.push_back(0.0);
	    }
          }
        }
      }
    }
  }

  cout << "write all Dirichlet faces\n";
  if(vDisp_x_Idx.size() > 0)
  {
    geomechfile << "GMNODE_BCDISPX\n";
    for(int i = 0; i < vDisp_x_Idx.size(); ++i)
    {
      geomechfile <<  vDisp_x_Idx[i] + 1 << "\t";
      geomechfile << vDisp_x_Val[i] << endl;
    }
    geomechfile << "/\n\n";
  }
  if(vDisp_y_Idx.size() > 0)
  {
    geomechfile << "GMNODE_BCDISPY\n";
    for(int i = 0; i < vDisp_y_Idx.size(); ++i)
    {
      geomechfile <<  vDisp_y_Idx[i] + 1 << "\t";
      geomechfile << vDisp_y_Val[i] << endl;
    }
    geomechfile << "/\n\n";
  }
  if(vDisp_z_Idx.size() > 0)
  {
    geomechfile << "GMNODE_BCDISPZ\n";
    for(int i = 0; i < vDisp_z_Idx.size(); ++i)
    {
      geomechfile <<  vDisp_z_Idx[i] + 1 << "\t";
      geomechfile << vDisp_z_Val[i] << endl;
    }
    geomechfile << "/\n\n";
  }
  geomechfile.close();

 outstring =   outputPath_ + "gm_fractures.txt";
 geomechfile.open(outstring.c_str());
 set<int>::iterator itsetint;

 int counter = 0;
 if (pSim->nInternalBoundaryFaces > 0)
 {
    int nFractures_ = 0;
    cout << "write all fractured faces\n";

    geomechfile << "GMFACE_FRACTURE_TO_FLOWCELL\n";
    for(itsetint = pSim->setIdenticalInternalMarker.begin(); itsetint != pSim->setIdenticalInternalMarker.end(); itsetint++, nFractures_++)
    {
      for (int i = 0; i < pSim->nPhysicalFacets; i++)
      {
        if( pSim->vsPhysicalFacet[i].nmark == *itsetint )
        {
	  geomechfile << pSim->vsPhysicalFacet[i].nface + 1 << "\t";
          geomechfile << pSim->vsPhysicalFacet[i].nfluid + 1 << endl;
          if( pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].nNeighbors !=2 )
          {
            cout << "Fracture interface # " << nFractures_ << endl;
            cout << "Global interface   # " << pSim->vsPhysicalFacet[i].nface << endl;
            cout << "Number od neighbors  " << pSim->vsFaceCustom[ pSim->vsPhysicalFacet[i].nface ].nNeighbors << endl;
            cout << "Wrong msh file. Mesh verticies are not connected on fracture interface" << endl;
            exit(0);
          }
        }
      }
    }
    geomechfile << "/" << endl << endl;

    nFractures_ = 0;
    geomechfile << "GMFACE_FRACTURE_CONDUCTIVITY" << endl;
    for(itsetint = pSim->setIdenticalInternalMarker.begin(); itsetint != pSim->setIdenticalInternalMarker.end(); itsetint++, nFractures_++)
    {
      for (int i = 0; i < pSim->nPhysicalFacets; i++)
      {
        if( pSim->vsPhysicalFacet[i].nmark == *itsetint )
	  geomechfile << pSim->vsFaceCustom[pSim->vsPhysicalFacet[i].nface].conductivity << endl;
      }
    }
    geomechfile << "/" << endl << endl;

    nFractures_ = 0;
    geomechfile << "GMFACE_FRACTURE_REGION\n";
    for(itsetint = pSim->setIdenticalInternalMarker.begin(); itsetint != pSim->setIdenticalInternalMarker.end(); itsetint++, nFractures_++)
    {
      for (int i = 0; i < pSim->nPhysicalFacets; i++)
      {
        if( pSim->vsPhysicalFacet[i].nmark == *itsetint )
        {
	  geomechfile << 1 << endl;
        }
      }
    }
    geomechfile << "/" << endl << endl;

    nFractures_ = 0;
    geomechfile << "GMFACE_FRACTURE_GROUP\n";
    for(itsetint = pSim->setIdenticalInternalMarker.begin(); itsetint != pSim->setIdenticalInternalMarker.end(); itsetint++, nFractures_++)
    {
      for (int i = 0; i < pSim->nPhysicalFacets; i++)
      {
        if( pSim->vsPhysicalFacet[i].nmark == *itsetint )
        {
	  geomechfile << nFractures_ + 1 << endl;
        }
      }
    }
    geomechfile << "/" << endl << endl;
    geomechfile.close();
 }
 geomechfile.close();

  outstring =   outputPath_ + "fl_pres.txt";
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

  outstring =   outputPath_ + "fl_temp.txt";
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

  outstring =   outputPath_ + "fl_zmf.txt";
 geomechfile.open(outstring.c_str());

  geomechfile << "ZMF\n";
  for ( int k = 0; k < 5; ++k )
  {
    for ( int iface = 0; iface < pSim->nFaces; iface++ )
    {
      // internal suface
      if ( pSim->vsFaceCustom[iface].nMarker > 0 )
      {
        int n1_ = pSim->vsFaceCustom[iface].vNeighbors[0];
        int n2_ = pSim->vsFaceCustom[iface].vNeighbors[1];
        geomechfile << min ( pSim->vsCellRockProps[n1_].zmf[k], pSim->vsCellRockProps[n2_].zmf[k] ) << endl;
      }
    }

    for(int ib = 0; ib < pSim->nCells; ++ib)
      geomechfile << pSim->vsCellRockProps[ib].zmf[k] << endl;
  }

  geomechfile << "\n/\n\n";
  geomechfile.close();


#if 0
  outstring =   outputPath_ + "gm_fracprops.txt";
  geomechfile.open ( outstring.c_str() );

  geomechfile <<  "GMCONTACT_THERMAL_EXPANSION" << endl;
  geomechfile <<  0.0 << endl;
  geomechfile << "/" << endl << endl;
  vector<double> a_;
  vector<double> b_;
  vector<double> c_;
  vector<double> d_;
  vector<double> e_;
  vector<double> f_;

  // Default Tengiz values
  a_.push_back(-5e-4); b_.push_back(0);
  a_.push_back(-5e-5); b_.push_back(0);
  a_.push_back(-5e-6); b_.push_back(0);
  a_.push_back(0.0); b_.push_back(0.0);
  for(unsigned int i = 4; i < 13; i++)
  {
    a_.push_back(0); b_.push_back(b_[i-1] + 25);
  }
  for(unsigned int i = 13; i < 20; i++)
  {
    a_.push_back(0); b_.push_back(b_[i-1] + 50);
  }

  c_.assign(b_.size(),0);
  for(unsigned int i = 0; i < b_.size(); i++)
    c_[i] = b_[i] * 0.6 + 100.0;

  d_.assign(b_.size(),0);
  e_.assign(b_.size(),0);
  f_.assign(b_.size(),0);
  for(unsigned int i = 3; i < 20; i++)
  {
   d_[i] = 7.5 / (1.0 + exp(1.5e-2*(b_[i]+10))) + 1;
   e_[i] = d_[i];
   f_[i] = 4.0 / (1.0 + exp(1.0e-2 * b_[i]))+0.81;
  }
  for(unsigned int i = 0; i < 3; i++)
  {
   d_[i] = d_[3] - 1e15 * a_[i] * a_[i] * a_[i] / 12./10.; //10 is a reference conductivity
   e_[i] = d_[i];
   f_[i] = f_[3] - a_[i] / 5e-3; //1e-3 is a refernce aperture
  }
  // Default Tengiz values

  geomechfile << "GMCONTACT_SLIP_RAMP" << endl;
  geomechfile << 5e-4 << "\t" << 1e-3 << "\t" << 1.0 << "\t" << 1.0 << endl;
  geomechfile << "/" << endl << endl;

  int conductivity_table = 1;
  if (conductivity_table < 0 )
  {
    for(unsigned int i = 3; i < 20; i++)
     d_[i] = 150 / (1.0 + exp(2.1e-2*(b_[i]+90))) + 1;
    for(unsigned int i = 0; i < 3; i++)
     d_[i] = d_[3] - 1e15 * a_[i] * a_[i] * a_[i] / 12./10.; //10 is a reference conductivity
    for(unsigned int i = 0; i < 20; i++)
     e_[i] = d_[i];
  }
  else if (conductivity_table == 0 )
  {
     d_.assign(d_.size(),1);
     e_.assign(d_.size(),1);
  }
  else
  {
    double c = double(conductivity_table);
    double a = (c - 0.9999*c) / (b_[3] - b_[b_.size()-1]);
    double b = c - a * b_[3];
    e_.assign(d_.size(),0);
    for(unsigned int i = 3; i < e_.size(); i++)
      e_[i] = (a * b_[i] + b) * d_[i];

    for(unsigned int i = 0; i < 3; i++)
    {
      d_[i] = d_[3] + fabs(a_[i] * a_[i] * a_[i]) / 12.0 / 1e-15 / 10.;
      e_[i] = e_[3] + fabs(a_[i] * a_[i] * a_[i]) / 12.0 / 1e-15 / 10.;
    }
  }

  int volume_table = 0;
  if (volume_table < 1 )
     f_.assign(d_.size(),1);

  double slip_coefficient = 0.6;
  geomechfile <<  "GMCONTACT_NORMAL_PROPS" << endl;
  for ( int i = 0; i < a_.size(); ++i )
  {
    geomechfile << a_[i] << "\t";
    geomechfile << b_[i] << "\t";
    geomechfile << b_[i] * slip_coefficient << "\t";
    geomechfile << d_[i] << "\t";
    geomechfile << e_[i] << "\t";
    geomechfile << f_[i] << endl;
  }
  geomechfile << "/" << endl << endl;
  geomechfile.close();
#endif

  outstring =   outputPath_ + "fl_satnum.txt";
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


 outstring =   outputPath_ + "fl_thc.txt";
    geomechfile.open(outstring.c_str());

  geomechfile << "THCROCK" << endl;
  for(int iface = 0; iface < pSim->nFaces; iface++)
  {
    // internal suface
    if( pSim->vsFaceCustom[iface].nMarker > 0)
    {
      int n1_ = pSim->vsFaceCustom[iface].vNeighbors[0];
      int n2_ = pSim->vsFaceCustom[iface].vNeighbors[1];
      geomechfile << min(pSim->vsCellRockProps[n1_].thc, pSim->vsCellRockProps[n2_].thc) << endl;
    }
  }
  for(int ib = 0; ib < pSim->nCells; ++ib)
    geomechfile << pSim->vsCellRockProps[ib].thc << endl;

  geomechfile << "/" << endl;
  geomechfile.close();


  cout << "write all wells\n";
  if ( pSim->nWells > 0 )
  {
    outstring =   outputPath_ + "fl_wells.txt";
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
