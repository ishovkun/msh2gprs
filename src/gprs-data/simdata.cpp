#include "simdata.hpp"


#include "Point.hpp"
// #include "GJK_Algorithm.hpp"
#include "Rectangle.hpp"
#include "CollisionGJK.hpp"
#include <muParser.h>

#define SPECIAL_CELL = 999
#include <algorithm>

const std::size_t EMBEDDED_FRACTURE_CELL = 9999992;

using Point = angem::Point<3, double>;

SimData::SimData(const string & inputstream, const SimdataConfig & config)
    :
    config(config)
{
  pStdElement = new StandardElements();

  instream = inputstream;

  outstream = inputstream;
  stringstream streamword(inputstream);
  streamword.imbue(locale(locale(), new wordtokens()));
  istream_iterator<string> begin(streamword);
  istream_iterator<string> end;
  vector<string> vwords(begin, end);
  outstream = vwords[0];
  nNodes = 0;
  nCells = 0;
//
  dNotNumber = -999.999;

  vPointPass.push_back(0);
  vPointPass.push_back(2);
  vPointPass.push_back(1);

  vPointCoord.resize(3, 0.0);

  vvPlate.resize(3);
  for (int i = 0; i < 3; i++) vvPlate[i].resize(3, 0.0);

  // Boundary conditions and initial state
  Sxx = 0.0;// 399.0= 0.748 * Szz bar [default: 0]
  Syy = 0.0; // 424.0 = 0.795 * Szz bar [default: 300.0]
  Szz = 0.0; // 533 = 2172*9.81*2500*1e-5 (rho*g*h) bar [default: 600]
  Syz = 0.0;
  Sxz = 0.0;
  Sxy = 0.0;

  //wells
  nWells = 2;
  vsWell.resize(nWells);

  // well 1 at the footwall block
  vsWell[0].vWellCoordinate.clear();
  vsWell[0].vWellCoordinate.push_back(-250.0); // x
  vsWell[0].vWellCoordinate.push_back(-0.0); // 45-30-90 model
  vsWell[0].vWellCoordinate.push_back(-1450.0); // z0
  vsWell[0].vWellCoordinate.push_back(-1550.0); // z1
  vsWell[0].Type = "WCONPROD";
  vsWell[0].radius_poisk = 6.0; // m

  // well 2 at the hangingwall block
  vsWell[1].vWellCoordinate.clear();
  vsWell[1].vWellCoordinate.push_back(250.0); // x
  vsWell[1].vWellCoordinate.push_back(-0.0); // 45-30-90 model
  vsWell[1].vWellCoordinate.push_back(-1450.0); // z0
  vsWell[1].vWellCoordinate.push_back(-1550.0); // z1
  vsWell[1].Type = "WCONPROD";
  vsWell[1].radius_poisk = 6.0; // m

  // Kirill's renumbering
  pRenum = new renum();
}

SimData::~SimData()
{
}


void SimData::defineEmbeddedFractureProperties()
{
  std::size_t ef_ind = 0;
  // class that checks if shapes collide
  angem::CollisionGJK<double> collision;
  // non-const since fracture is adjusted to avoid collision with vertices
  for (auto & frac_conf : config.fractures)
  {
    vsEmbeddedFractures.emplace_back();
    auto & frac = vsEmbeddedFractures.back();
    Point total_shift = {0, 0, 0};

 redo_collision:
    // find cells intersected by the fracture
    frac.cells.clear();
    for(std::size_t icell = 0; icell < nCells; icell++)
    {
      const auto & cell = vsCellCustom[icell];

      std::vector<angem::Point<3,double> *> verts;
      for (const auto & ivertex : cell.vVertices)
        verts.push_back(&vvVrtxCoords[ivertex]);

      angem::Shape<double> pcell(verts);

      if (collision.check(frac_conf.body, pcell))
      {
        frac.cells.push_back(icell);

        // check if some vertices are too close to the fracture
        // and if so move a fracture a little bit
        for (const auto & ivertex : cell.vVertices)
        {
          const auto & vert = vvVrtxCoords[ivertex];
          const auto vc = vsCellCustom[icell].center - vert;
          if ( fabs(frac_conf.body.plane.distance(vert)/vc.norm()) < 1e-6 )
          {
            const auto shift = 1e-5 * frac_conf.body.plane.normal();
            total_shift += shift;
            std::cout << "shifting fracture: " << shift ;
            std::cout << " due to collision with vertex: " << ivertex;
            std::cout << std::endl;
            // std::cout << frac.plane.point << "\t";
            frac_conf.body.move(shift);
            goto redo_collision;
          }
        }  // end adjusting fracture
      }  // end collision processing
    }    // end cell loop
    std::cout << "Total shift = " << total_shift << std::endl;
    const std::size_t n_efrac_cells = frac.cells.size();
    if (n_efrac_cells == 0)
    {
      vsEmbeddedFractures.pop_back();
      continue;
    }

    // fill out properties
    vsEmbeddedFractures[ef_ind].points.assign(n_efrac_cells,
                                              frac_conf.body.center());
    vsEmbeddedFractures[ef_ind].dip.assign(n_efrac_cells,
                                           frac_conf.body.plane.dip_angle());
    vsEmbeddedFractures[ef_ind].strike.assign(n_efrac_cells,
                                              frac_conf.body.plane.strike_angle());

    vsEmbeddedFractures[ef_ind].cohesion = frac_conf.cohesion;
    vsEmbeddedFractures[ef_ind].friction_angle = frac_conf.friction_angle;
    vsEmbeddedFractures[ef_ind].dilation_angle = frac_conf.dilation_angle;

    ef_ind++;
  }  // end efracs loop

}


std::size_t SimData::n_default_vars()
{
  SimdataConfig dummy;
  // return 0;
  return dummy.all_vars.size();
}


void SimData::defineRockProperties()
{
  // print header
  std::cout << "function parsers setup" << std::endl;
  std::cout << "Variables:" << std::endl;
  for (std::size_t i=0; i<config.all_vars.size(); ++i)
  {
    std::cout << config.all_vars[i] << "\t";
    if ((i+1)%10 == 0)
      std::cout << std::endl;
  }
  std::cout << std::endl;

  // get number of variables in an empty simdataconfig - should be 3=x+y+z
  const std::size_t shift = n_default_vars();

  // resize rock properties
  vsCellRockProps.resize(nCells);

  const std::size_t n_variables = config.all_vars.size();
  std::vector<double> vars(n_variables);

  // loop various domain configs:
  // they may have different number of variables and expressions
  for (const auto & conf: config.domains)
  {
    const std::size_t n_expressions = conf.expressions.size();
    std::vector<mu::Parser> parsers(n_expressions);

    // assign expressions and variables
    for (std::size_t i=0; i<n_expressions; ++i)
    {
      for (std::size_t j=0; j<n_variables; ++j)
      {
        try {
          parsers[i].DefineVar(config.all_vars[j], &vars[j]);
        }
        catch(mu::Parser::exception_type & e)
        {
          std::cout << _T("Initialization error:  ") << e.GetMsg() << endl;
          std::cout << "when setting variable '"
                    << config.all_vars[j]
                    << "'" << std::endl;
        }

      }
      try {
        parsers[i].SetExpr(conf.expressions[i]);
      }
      catch(mu::Parser::exception_type & e)
      {
        std::cout << _T("Initialization error:  ") << e.GetMsg() << endl;
        std::cout << "when setting expression '"
                  << conf.expressions[i]
                  << "'" << std::endl;
      }
    }

    // loop cells and evaluate expressions
    for ( int icell = 0; icell < nCells; icell++ )
      if ( vsCellCustom[icell].nMarker == conf.label ) // Regular cells
      {
        std::fill(vars.begin(), vars.end(), 0);
        vars[0] = vsCellCustom[icell].center(0);  // x
        vars[1] = vsCellCustom[icell].center(1);  // y
        vars[2] = vsCellCustom[icell].center(2);  // z

        // Evaluate expression -> write into variable
        for (std::size_t i=0; i<n_expressions; ++i)
          vars[conf.local_to_global_vars.at(i)] = parsers[i].Eval();

        // copy vars to cell properties
        vsCellRockProps[icell].v_props.resize(n_variables-shift);
        // start from 3 to skip x,y,z
        for (std::size_t j=shift; j<n_variables; ++j)
        {
          try
          {
            vsCellRockProps[icell].v_props[j-shift] = vars[j];
          }
          catch (std::out_of_range & e)
          {
            vsCellRockProps[icell].v_props[j-shift] = 0;
          }
        }
      }  // cell loop

  }  // end domain loop

}


void SimData::readGmshFile()
{
  Gelement Element3D, Element2D;

  fstream inputfile;
  inputfile.open (instream.c_str(), fstream::in);

  std::string inputline;
  getline (inputfile, inputline);

  std::stringstream streamline(stringstream::in | stringstream::out);
  streamline << inputline;
  streamline.imbue(std::locale(std::locale(), new tokens()));
  std::istream_iterator<std::string> begin;
  std::istream_iterator<std::string> end;
  std::vector<std::string> vstrings;

  copy( istream_iterator<string>(streamline), istream_iterator<string>(),back_inserter(vstrings) );

    while(vstrings[0] != "$Nodes")
    {
      getline(inputfile, inputline);
      streamline.clear(); vstrings.clear();
      streamline << inputline;
      streamline.imbue(std::locale(std::locale(), new tokens()));
      copy( istream_iterator<string>(streamline), istream_iterator<string>(),back_inserter(vstrings) );
    }

    getline(inputfile, inputline);
    streamline.clear(); vstrings.clear();
    streamline << inputline;
    streamline.imbue(std::locale(std::locale(), new tokens()));
    copy( istream_iterator<string>(streamline), istream_iterator<string>(),back_inserter(vstrings) );

    nNodes = atoi(vstrings[0].c_str());

    vvVrtxCoords.resize(nNodes);
    for(int i = 0; i < nNodes; i++ )
    {
      getline(inputfile, inputline);
      streamline.clear(); vstrings.clear();
      streamline << inputline;
      streamline.imbue(std::locale(std::locale(), new tokens()));
      copy( istream_iterator<string>(streamline), istream_iterator<string>(),back_inserter(vstrings) );

      for (int j = 0; j < 3; j++)
        vvVrtxCoords[i][j] = atof(vstrings[j+1].c_str());
    }
    nNodes = vvVrtxCoords.size();

    maxVrtxCoordsX = -1e10;
    maxVrtxCoordsY = -1e10;
    maxVrtxCoordsZ = -1e10;
    minVrtxCoordsX = 1e10;
    minVrtxCoordsY = 1e10;
    minVrtxCoordsZ = 1e10;

    for(int i = 0; i < nNodes; i++ )
    {
      if(vvVrtxCoords[i][0] > maxVrtxCoordsX ) maxVrtxCoordsX = vvVrtxCoords[i][0];
      if(vvVrtxCoords[i][1] > maxVrtxCoordsY ) maxVrtxCoordsY = vvVrtxCoords[i][1];
      if(vvVrtxCoords[i][2] > maxVrtxCoordsZ ) maxVrtxCoordsZ = vvVrtxCoords[i][2];

      if(vvVrtxCoords[i][0] < minVrtxCoordsX ) minVrtxCoordsX = vvVrtxCoords[i][0];
      if(vvVrtxCoords[i][1] < minVrtxCoordsY ) minVrtxCoordsY = vvVrtxCoords[i][1];
      if(vvVrtxCoords[i][2] < minVrtxCoordsZ ) minVrtxCoordsZ = vvVrtxCoords[i][2];
    }

    while(vstrings[0] != "$Elements")
    {
     getline(inputfile, inputline);
     streamline.clear(); vstrings.clear();
     streamline << inputline;
     streamline.imbue(std::locale(std::locale(), new tokens()));
     copy( istream_iterator<string>(streamline), istream_iterator<string>(),back_inserter(vstrings) );
    }

    getline(inputfile, inputline);
    streamline.clear(); vstrings.clear();
    streamline << inputline;
    streamline.imbue(std::locale(std::locale(), new tokens()));
    copy( istream_iterator<string>(streamline), istream_iterator<string>(),back_inserter(vstrings) );

    int n = atoi(vstrings[0].c_str());

    nCells = 0; nFaces = 0;
    for(int i = 0; i < n; i++ )
    {
      getline(inputfile, inputline);
      streamline.clear(); vstrings.clear();
      streamline << inputline;
      streamline.imbue(std::locale(std::locale(), new tokens()));
      copy( istream_iterator<string>(streamline), istream_iterator<string>(),back_inserter(vstrings) );

      int elem_type = atoi(vstrings[1].c_str());
      int ntags = atoi(vstrings[2].c_str());
      int node_list_run = 0;
      int node_list_end = 0;

      Element3D.nNeighbors = 0;
      Element3D.nVertices = 0;
      Element3D.center.clear();
      Element3D.vNeighbors.clear();
      Element3D.normal.clear();
      Element3D.vVertices.clear();
      Element3D.vVerticesNewnum.clear();
      Element3D.vVerticesSorted.clear();

      Element2D.nNeighbors = 0;
      Element2D.nVertices = 0;
      Element2D.center.clear();
      Element2D.vNeighbors.clear();
      Element2D.normal.clear();
      Element2D.vVertices.clear();
      Element2D.vVerticesNewnum.clear();
      Element2D.vVerticesSorted.clear();

        switch (elem_type)
        {
        case 2:
        node_list_run = 3 + ntags;
        node_list_end = node_list_run + 3;
            for (int j = node_list_run; j < node_list_end; j++) Element2D.vVertices.push_back(atoi(vstrings[j].c_str()));
            Element2D.vtkIndex = 5;
            Element2D.formIndex = TRGLE3;
            Element2D.nMarker = atoi(vstrings[3].c_str());
            Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
            vsFaceCustom.push_back( Element2D );
            nFaces++;
            break;

        case 3:
        node_list_run = 3 + ntags;
        node_list_end = node_list_run + 4;
            for (int j = node_list_run; j < node_list_end; j++) Element2D.vVertices.push_back(atoi(vstrings[j].c_str()));
            Element2D.vtkIndex = 9;
            Element2D.formIndex = QUAD4;
            Element2D.nMarker = atoi(vstrings[3].c_str());
            Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
            vsFaceCustom.push_back( Element2D );
            nFaces++;
            break;

        case 4:
        node_list_run = 3 + ntags;
        node_list_end = node_list_run + 4;
            for (int j = node_list_run; j < node_list_end; j++) Element3D.vVertices.push_back(atoi(vstrings[j].c_str()));
            Element3D.vtkIndex = 10;
            Element3D.formIndex = TETRA4;
            Element3D.nMarker = atoi(vstrings[3].c_str());
            vsCellCustom.push_back( Element3D );
            nCells++;
            break;

        case 5:
        node_list_run = 3 + ntags;
        node_list_end = node_list_run + 8;
            for (int j = node_list_run; j < node_list_end; j++) Element3D.vVertices.push_back(atoi(vstrings[j].c_str()));
            Element3D.vtkIndex = 12;
            Element3D.formIndex = PRISM8;
            Element3D.nMarker = atoi(vstrings[3].c_str());
            vsCellCustom.push_back( Element3D );
            nCells++;
            break;

        case 6:
        node_list_run = 3 + ntags;
        node_list_end = node_list_run + 6;
            for (int j = node_list_run; j < node_list_end; j++) Element3D.vVertices.push_back(atoi(vstrings[j].c_str()));
            Element3D.vtkIndex = 13;
            Element3D.formIndex = PRISM6;
            Element3D.nMarker = atoi(vstrings[3].c_str());
            vsCellCustom.push_back( Element3D );
            nCells++;
            break;

        case 9:
        node_list_run = 3 + ntags;
        node_list_end = node_list_run + 6;
            for (int j = node_list_run; j < node_list_end; j++) Element2D.vVertices.push_back(atoi(vstrings[j].c_str()));
            Element2D.vtkIndex = 22;
            Element2D.formIndex = TRGLE6;
            Element2D.nMarker = atoi(vstrings[3].c_str());
            Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
            vsFaceCustom.push_back( Element2D );
            nFaces++;
            break;

        case 11:
        node_list_run = 3 + ntags;
        node_list_end = node_list_run + 10;
            for (int j = node_list_run; j < node_list_end; j++) Element3D.vVertices.push_back(atoi(vstrings[j].c_str()));

            // because different numbering (vtk vs gmsh)
            std::swap(Element3D.vVertices[8], Element3D.vVertices[9]);

            Element3D.vtkIndex = 24;
            Element3D.formIndex = TETRA10;
            Element3D.nMarker = atoi(vstrings[3].c_str());
            vsCellCustom.push_back( Element3D );
            nCells++;
            break;

        case 16:
        node_list_run = 3 + ntags;
        node_list_end = node_list_run + 8;
            for (int j = node_list_run; j < node_list_end; j++) Element2D.vVertices.push_back(atoi(vstrings[j].c_str()));
            Element2D.vtkIndex = 23;
            Element2D.formIndex = QUAD8;
            Element2D.nMarker = atoi(vstrings[3].c_str());
            Element2D.nMarker *= checkReservedBoundaryName(Element2D.nMarker);
            vsFaceCustom.push_back( Element2D );
            nFaces++;
            break;

        case 17:
        node_list_run = 3 + ntags;
        node_list_end = node_list_run + 20;
            for (int j = node_list_run; j < node_list_end; j++) Element3D.vVertices.push_back(atoi(vstrings[j].c_str()));
            Element3D.vtkIndex = 25;
            Element3D.formIndex = PRISM20;
            Element3D.nMarker = atoi(vstrings[3].c_str());
            vsCellCustom.push_back( Element3D );
            nCells++;
            break;

        case 18:
        node_list_run = 3 + ntags;
        node_list_end = node_list_run + 15;
            for (int j = node_list_run; j < node_list_end; j++) Element3D.vVertices.push_back(atoi(vstrings[j].c_str()));
            Element3D.vtkIndex = 26;
            Element3D.formIndex = PRISM15;
            Element3D.nMarker = atoi(vstrings[3].c_str());
            vsCellCustom.push_back( Element3D );
            nCells++;
            break;

        default:
            cout << "Element type " <<  elem_type << endl;
            cout << "Wrong element type. Supported: {2, 3, 4, 5, 6, 9, 11, 17, 18}\n";
            exit(-1);
            break;
        }
    }
    // c++ rule
    for(int icell = 0; icell < nCells; icell++)
    {
      vsCellCustom[icell].nVertices = vsCellCustom[icell].vVertices.size();
      for(int ivrtx = 0; ivrtx < vsCellCustom[icell].nVertices; ivrtx++ )
      {
        vsCellCustom[icell].vVertices[ivrtx] = vsCellCustom[icell].vVertices[ivrtx] - 1;
      }
    }

    for(int iface = 0; iface < nFaces; iface++)
    {
      vsFaceCustom[iface].nVertices = vsFaceCustom[iface].vVertices.size();
      for(int ivrtx = 0; ivrtx < vsFaceCustom[iface].nVertices; ivrtx++ )
      {
        vsFaceCustom[iface].vVertices[ivrtx] = vsFaceCustom[iface].vVertices[ivrtx] - 1;
      }
    }
    cout << "GMSH vrtxs " << nNodes << endl;
    cout << "GMSH cells " << nCells << endl;
    cout << "GMSH faces " << nFaces << endl;
  inputfile.close();
}

void SimData::extractInternalFaces()
{
  vector<int> vLocalPolygonVertices;
  set<int> setLocalPolygonVertices;
  stringstream vertices_stream;

  set<string>  setIdenticalPolygons;
  pair<set<string>::iterator, bool> pair_itstring_bool;

  // At first we should write all input polygons
  for(int iface = 0; iface < nFaces; iface++ )
  {
    setLocalPolygonVertices.clear();
    for(int ivrtx = 0; ivrtx < vsFaceCustom[iface].nVertices; ivrtx++ )
    {
      setLocalPolygonVertices.insert ( vsFaceCustom[iface].vVertices[ivrtx] );
    }
    set<int>::iterator it_set;
    vertices_stream.str ( "" );
    for(it_set = setLocalPolygonVertices.begin(); it_set != setLocalPolygonVertices.end(); ++it_set) vertices_stream << *it_set;
    //try write into set
    pair_itstring_bool = setIdenticalPolygons.insert ( vertices_stream.str() );
  }

  // we dont know how many polygons
  for ( int icell = 0; icell < nCells; icell++ )
  {
    int job_percent = int ( ( 100. * icell ) / ( nCells ) );
    cout << "\r    " << job_percent << "%";

    // loop by all element faces
    for ( int iface = 0; iface < pStdElement->elementProps[ vsCellCustom[icell].formIndex ].facesInElement; iface++ )
    {
      vLocalPolygonVertices.clear();
      setLocalPolygonVertices.clear();
      //num of nodes on face
      int nnodes_ = pStdElement->elementProps[ vsCellCustom[icell].formIndex].vvFacesNodes[iface].size();

      //write nodes into vector
      for ( int inode = 0; inode < nnodes_; inode++ )
      {
        int local_node = pStdElement->elementProps[ vsCellCustom[icell].formIndex].vvFacesNodes[iface][inode];
        int global_node = vsCellCustom[icell].vVertices[local_node];

        vLocalPolygonVertices.push_back ( global_node );
        setLocalPolygonVertices.insert ( global_node );
      }

      set<int>::iterator it_set;
      vertices_stream.str ( "" );

      for ( it_set = setLocalPolygonVertices.begin(); it_set != setLocalPolygonVertices.end(); ++it_set )
        vertices_stream << *it_set;

      //try write into set
      pair_itstring_bool = setIdenticalPolygons.insert ( vertices_stream.str() );

      //element will be written if we see it first time
      Gelement temporaryElement;
      if ( pair_itstring_bool.second == true )
      {
        // extend vFaceCustom
        temporaryElement.fluidElement = -1;
        temporaryElement.nNeighbors = 0;
        temporaryElement.nVertices = vLocalPolygonVertices.size();
        temporaryElement.vVertices.resize ( temporaryElement.nVertices );
        temporaryElement.vVertices = vLocalPolygonVertices;

        if ( temporaryElement.nVertices == 3 )
        {
          temporaryElement.vtkIndex = 5;
          temporaryElement.formIndex = TRGLE3;
        }

        if ( temporaryElement.nVertices == 6 )
        {
          temporaryElement.vtkIndex = 22;
          temporaryElement.formIndex = TRGLE6;
        }

        if ( temporaryElement.nVertices == 4 )
        {
          temporaryElement.vtkIndex = 9;
          temporaryElement.formIndex = QUAD4;
        }

        if ( temporaryElement.nVertices == 8 )
        {
          temporaryElement.vtkIndex = 23;
          temporaryElement.formIndex = QUAD8;
        }

        temporaryElement.nMarker = 0;

        vsFaceCustom.push_back ( temporaryElement );
      }
    }
  }
  nFaces = vsFaceCustom.size();

}

void SimData::convertGmsh2Sim()
{
  int counter_ = 0;
  for(int iface = 0; iface < nFaces; iface++)
  {
    methodElementCenter(iface, vsFaceCustom);
    if(vsFaceCustom[iface].nMarker != 0) methodFaceNormalVector(iface, vsFaceCustom);
  }
  methodChangeFacesNormalVector();

  counter_ = 0;
  corner_cell = 0;
  double minx_ = 1e-16;
  double miny_ = 1e-16;
  for(int icell = 0; icell < nCells; icell++)
  {
    methodElementCenter(icell, vsCellCustom);
    if( vsCellCustom[icell].center[0] < minx_ &&  vsCellCustom[icell].center[1] < miny_)
    {
       minx_ = vsCellCustom[icell].center[0];
       miny_ = vsCellCustom[icell].center[1];
       corner_cell = icell;
    }
  }

  cout << "\t sort verticies" << endl;
  for(int iface = 0; iface < nFaces; iface++)
  {
    vsFaceCustom[iface].vVerticesSorted.resize( vsFaceCustom[iface].nVertices );
    vsFaceCustom[iface].vVerticesSorted = vsFaceCustom[iface].vVertices;
    sort(vsFaceCustom[iface].vVerticesSorted.begin(), vsFaceCustom[iface].vVerticesSorted.end());
  }

  for(int icell = 0; icell < nCells; icell++)
  {
    vsCellCustom[icell].vVerticesSorted.resize( vsCellCustom[icell].nVertices );
    vsCellCustom[icell].vVerticesSorted = vsCellCustom[icell].vVertices;
    sort(vsCellCustom[icell].vVerticesSorted.begin(), vsCellCustom[icell].vVerticesSorted.end());
  }

  cout << "\t find face neighbor cells / cell neighbor faces (slow)" << endl;
  vsetPolyhedronPolygon.resize(nCells);
  vsetPolygonPolyhedron.resize(nFaces);
  for ( int iface = 0; iface < nFaces; iface++ )
  {
    int job_percent = int ( ( 100. * iface ) / ( nFaces ) );
    cout << "\r    " << job_percent << "%";
    for ( int icell = 0; icell < nCells; icell++ )
    {
      if ( includes ( vsCellCustom[icell].vVerticesSorted.begin(), vsCellCustom[icell].vVerticesSorted.end(),
                      vsFaceCustom[iface].vVerticesSorted.begin(), vsFaceCustom[iface].vVerticesSorted.end() ) )
      {
        vsetPolyhedronPolygon[icell].insert ( iface );
        vsetPolygonPolyhedron[iface].insert ( icell );
      }
    }
  }

  for ( int iface = 0; iface < nFaces; iface++ )
  {
    if(vsetPolygonPolyhedron[iface].size() == 0)
    {
      cout << endl << "Polygon " << iface << " (" << nFaces << ") has no connected polyhedrons!" << endl;
      exit(0);
    }
  }

  for ( int ic = 0; ic < nCells; ic++ )
  {
    if(vsetPolyhedronPolygon[ic].size() == 0)
    {
      cout << endl << "Polyhedron " << ic << " has no connected polygons!" << endl;
      exit(0);
    }
  }

  cout << "\t find cell neighbor cells " << endl;
  vector<set<int> > vsetCellNeighborCells;
  vsetCellNeighborCells.resize(nCells);
  for ( int icell = 0; icell < nCells; icell++ )
  {
    set<int>::iterator it_set, it_set_f;

    for ( it_set = vsetPolyhedronPolygon[icell].begin(); it_set != vsetPolyhedronPolygon[icell].end(); ++it_set )
    {
      int nface = *it_set;

      for ( it_set_f = vsetPolygonPolyhedron[nface].begin(); it_set_f != vsetPolygonPolyhedron[nface].end(); ++it_set_f )
      {
        int ncell = *it_set_f;
        if(ncell != icell) vsetCellNeighborCells[icell].insert(ncell);
      }
    }

  }

  for ( int icell = 0; icell < nCells; icell++ )
  {
    vsCellCustom[icell].vNeighbors.clear();
    set<int>::iterator it_set;
    for ( it_set = vsetCellNeighborCells[icell].begin(); it_set != vsetCellNeighborCells[icell].end(); ++it_set )
    {
      vsCellCustom[icell].vNeighbors.push_back( *it_set );
    }

    vsCellCustom[icell].nNeighbors = vsCellCustom[icell].vNeighbors.size();
  }


  counter_ = 0;
  for ( int icell = 0; icell < nCells; icell++ )
  {

    for ( int i = 0; i < vsCellCustom[icell].nVertices; i++ )
    {
      vsCellCustom[icell].vVerticesNewnum.push_back(counter_);
      counter_++;
    }
  }

  cout << "\t choise support cell for internal facets " << endl;
  vector<double> vDatumDistance;
  for ( int iface = 0; iface < nFaces; iface++ )
  {
    vsFaceCustom[iface].vNeighbors.clear();

    set<int>::iterator it_set;

    for ( it_set = vsetPolygonPolyhedron[iface].begin(); it_set != vsetPolygonPolyhedron[iface].end(); ++it_set )
      vsFaceCustom[iface].vNeighbors.push_back ( *it_set );

    vsFaceCustom[iface].nNeighbors = vsFaceCustom[iface].vNeighbors.size();

    if( vsFaceCustom[iface].nNeighbors == 0 )
    {
      cout << "Face " << iface << " has no connected cells!" << endl;
      exit(0);
    }

    /// choise support cell for internal facets
    if ( vsFaceCustom[iface].nNeighbors == 2 && vsFaceCustom[iface].nMarker > 0 )
    {
      double cosa = 0;
      for ( int idx = 0; idx < 3; idx++ )
        cosa += (vsCellCustom[ vsFaceCustom[iface].vNeighbors[0] ].center[idx] - vsFaceCustom[iface].center[idx]) * vsFaceCustom[iface].normal[idx];

      if(cosa > 0) swap(vsFaceCustom[iface].vNeighbors[0], vsFaceCustom[iface].vNeighbors[1]);
    }

  }

}

void SimData::methodElementCenter(int nelem, vector<Gelement> &vsElement)
{
  int nodes_in_elem = vsElement[nelem].nVertices;

  vsElement[nelem].center.clear();

  for (int inodes = 0; inodes < nodes_in_elem; inodes++)
  {
    int gl_node_num = vsElement[nelem].vVertices[inodes];
    for (int idx = 0; idx < 3; idx++)
    {
      vsElement[nelem].center[idx] += vvVrtxCoords[gl_node_num][idx] / nodes_in_elem;
    }
  }

  vsElement[nelem].center_distance = 0.0;
  for (int idx = 0; idx < 3; idx++)
  {
    vsElement[nelem].center_distance += vsElement[nelem].center[idx] * vsElement[nelem].center[idx];
  }
  vsElement[nelem].center_distance = sqrt(vsElement[nelem].center_distance);

  double distance, buf;
  int node_master, node_slave;
  switch ( vsElement[nelem].formIndex )
  {
  case TETRA4:
    distance = 0.0;
    node_master = vsElement[nelem].vVertices[0];
    for ( int inode = 1; inode < nodes_in_elem; ++inode )
    {
      node_slave = vsElement[nelem].vVertices[inode];
      buf = 0.0;
      for ( int idx = 0; idx < 3; idx++ )
      {
        buf += ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] ) * ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] );
      }
      distance += sqrt ( buf );

    }
    vsElement[nelem].thickness = distance / (nodes_in_elem-1);
    break;
  case TETRA10:
    distance = 0.0;
    node_master = vsElement[nelem].vVertices[0];
    for ( int inode = 1; inode < nodes_in_elem; ++inode )
    {
      node_slave = vsElement[nelem].vVertices[inode];
      buf = 0.0;
      for ( int idx = 0; idx < 3; idx++ )
      {
        buf += ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] ) * ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] );
      }
      distance += sqrt ( buf );

    }
    vsElement[nelem].thickness = distance / (nodes_in_elem-1);
    break;
  case PRISM6:
    distance = 0.0;
    for ( int inode = 0; inode < nodes_in_elem / 2; ++inode )
    {
      node_master = vsElement[nelem].vVertices[inode];
      node_slave = vsElement[nelem].vVertices[inode + nodes_in_elem / 2];
      buf = 0.0;
      for ( int idx = 0; idx < 3; idx++ )
      {
        buf += ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] ) * ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] );
      }
      distance += sqrt ( buf );
    }
    vsElement[nelem].thickness = distance / (nodes_in_elem / 2);
    break;
  case PRISM8:
    distance = 0.0;
    for ( int inode = 0; inode < nodes_in_elem / 2; ++inode )
    {
      node_master = vsElement[nelem].vVertices[inode];
      node_slave = vsElement[nelem].vVertices[inode + nodes_in_elem / 2];
      buf = 0.0;
      for ( int idx = 0; idx < 3; idx++ )
      {
        buf += ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] ) * ( vvVrtxCoords[node_master][idx] - vvVrtxCoords[node_slave][idx] );
      }
      distance += sqrt ( buf );
    }
    vsElement[nelem].thickness = distance / (nodes_in_elem / 2);
    break;
  default:
    break;
  }
}


void SimData::methodFaceNormalVector(int nelem, vector<Gelement> &vsElement)
{

  for (int i = 0; i < 3; i++) vvPlate[i].assign(3, 0.0);

  for (int inodes = 0; inodes < 3; inodes++)
  {
    int gl_node_num = vsElement[nelem].vVertices[vPointPass[inodes]];
    for (int idx = 0; idx < 3; idx++)
    {
      vvPlate[inodes][idx] = vvVrtxCoords[gl_node_num][idx];
    }
  }
  // calculate basis
  vPointCoord[0] = vvPlate[1][1] * vvPlate[2][2] - vvPlate[1][2] * vvPlate[2][1] +
                   vvPlate[2][1] * vvPlate[0][2] - vvPlate[0][1] * vvPlate[2][2] +
                   vvPlate[0][1] * vvPlate[1][2] - vvPlate[1][1] * vvPlate[0][2];

  vPointCoord[1] = vvPlate[0][0] * vvPlate[2][2] - vvPlate[0][0] * vvPlate[1][2] +
                   vvPlate[1][0] * vvPlate[0][2] - vvPlate[1][0] * vvPlate[2][2] +
                   vvPlate[2][0] * vvPlate[1][2] - vvPlate[2][0] * vvPlate[0][2];

  vPointCoord[2] = -vvPlate[1][0] * vvPlate[0][1] + vvPlate[0][0] * vvPlate[1][1] -
                   vvPlate[0][0] * vvPlate[2][1] + vvPlate[1][0] * vvPlate[2][1] +
                   vvPlate[2][0] * vvPlate[0][1] - vvPlate[2][0] * vvPlate[1][1];

  //calculate length
  double length = 0.0;
  for (int i = 0; i < 3; i++) length += vPointCoord[i] * vPointCoord[i];

  length = sqrt(length);

  // calculate normal basis
  vsElement[nelem].normal.clear();

  for (int i = 0; i < 3; i++) vsElement[nelem].normal[i] = vPointCoord[i] / length;


}

void SimData::methodChangeFacesNormalVector()
{
  pair<set<int>::iterator, bool> pairIterBool;
  angem::Point<3,double> vDatumNormal = {0,0,0};

  vector<int> vFacevVertices;

  nExternalBoundaryFaces = 0;
  nInternalBoundaryFaces = 0;

  for(int iface = 0; iface < nFaces; iface++)
  {
    if( vsFaceCustom[iface].nMarker < 0)
    {
      setIdenticalExternalMarker.insert( vsFaceCustom[iface].nMarker );
      nExternalBoundaryFaces++;
    }

    if( vsFaceCustom[iface].nMarker > 0)
    {
      pairIterBool = setIdenticalInternalMarker.insert( vsFaceCustom[iface].nMarker );
      double a = 1e-3;
      vIdenticalInternalFacetPerm.push_back ( a*a/12/1e-15 );
      vIdenticalInternalFacetAperture.push_back ( 0.001 );
      vIdenticalInternalFacetFFpermMult.push_back ( 1.0 );

      nInternalBoundaryFaces++;
    }
  }
  set<int>::iterator itintset;
  for (itintset = setIdenticalInternalMarker.begin(); itintset != setIdenticalInternalMarker.end(); ++itintset)
  {
    // datum normal vector
    for(int iface = 0; iface < nFaces; iface++)
    {
      if( vsFaceCustom[iface].nMarker == *itintset)
      {
        vDatumNormal = vsFaceCustom[iface].normal;
        break;
      }
    }

    for(int iface = 0; iface < nFaces; iface++)
    {
      if( vsFaceCustom[iface].nMarker == *itintset)
      {
        double cosa = 0;
        for(int idx = 0; idx < 3; idx++) cosa += vDatumNormal[idx] * vsFaceCustom[iface].normal[idx];

        // non collinear vector. change verticies order
        if(cosa < 0.0)
        {
          vFacevVertices.clear();
          for(int ivrtx = 0; ivrtx < vsFaceCustom[iface].nVertices; ivrtx++)
          {
            vFacevVertices.push_back( vsFaceCustom[iface].vVertices[ vsFaceCustom[iface].nVertices - ivrtx - 1] );
          }
          vsFaceCustom[iface].vVertices.swap(vFacevVertices);

          for(int idx = 0; idx < 3; idx++) vsFaceCustom[iface].normal[idx] *= -1.0;

        }

      }
    }
  }

}

void SimData::methodRandomRockProperties()
{
 for (int i = 0; i < nCells; i++)
 {
   double x = vsCellCustom[i].center[0];
   double y = vsCellCustom[i].center[1];

   double lognormal_value = createLognormalDistribution(0.01, 0.001);
   double value = lognormal_value *
                  1.0/(1.0+exp((x-200) / 5)) * 1.0/(1.0+exp(-(x-20) / 5)) *
                  1.0/(1.0+exp((y-200) / 5)) * 1.0/(1.0+exp(-(y-20) / 5));

   vsCellRockProps[i].perm = value;

   lognormal_value = createLognormalDistribution(2.0, 0.0001);
   value = lognormal_value;
   vsCellRockProps[i].young = value;
 }
}

double SimData::createLognormalDistribution(double E, double S)
{
  vector<double> vVars;
  double value = 0.0;

  for(int i = 0; i < 12; i++)
  {
    vVars.push_back( (rand() % 100 + 1) / 100.);
  }

  for(int i = 0; i < 12; i++)
  {
    value += vVars[i];
  }
  value -= 6.0;

  return(E + S * value);
}


void SimData::splitInternalFaces()
{
  // int counter_;
  // cout << endl << "Create set of stick vertices (slow)" << endl;
  // vector<set<int> > vsetGlueVerticies;
  // counter_ = 0;
  // for ( int icell = 0; icell < nCells; icell++ )
  // {
  //   counter_ += vsCellCustom[icell].nVertices;
  // }
  // vsetGlueVerticies.resize(counter_);
  // for(int i = 0; i < counter_; i++) vsetGlueVerticies[i].insert(i);

  // // TODO
  // vector<bool> exlude_some_verticies(nNodes,false);
  // for ( int iface = 0; iface < nFaces; iface++ )
  // {
  //     vector<int>::iterator it_face_vrtx;
  //     for ( it_face_vrtx = vsFaceCustom[ iface ].vVertices.begin();  it_face_vrtx != vsFaceCustom[ iface ].vVertices.end(); ++it_face_vrtx)
  //     {
  //       if(vsFaceCustom[ iface ].nMarker == -3333332 || vsFaceCustom[ iface ].nMarker == -3333331)
  //         exlude_some_verticies[*it_face_vrtx] = true;
  //     }
  // }

  // vector<double> vVerticesPair; vVerticesPair.resize(2, 0);
  // for ( int iface = 0; iface < nFaces; iface++ )
  // {
  //   int job_percent = int ( ( 100. * iface ) / ( nFaces ) );
  //   cout << "\r    " << job_percent << "%";
  //   /// non phisical face
  //   if ( vsFaceCustom[ iface ].nMarker == 0 && vsFaceCustom[ iface ].nNeighbors == 2 )
  //   {
  //     vector<int>::iterator it_face_vrtx;

  //     /// loop by all face vertices
  //     for ( it_face_vrtx = vsFaceCustom[ iface ].vVertices.begin();  it_face_vrtx != vsFaceCustom[ iface ].vVertices.end(); ++it_face_vrtx)
  //     {
  //       /// loop by neighbor cells
  //       int n_polyhedron = 0;
  //       set<int>::iterator it_polyhedron;
  //       for ( it_polyhedron = vsetPolygonPolyhedron[iface].begin(); it_polyhedron != vsetPolygonPolyhedron[iface].end(); ++it_polyhedron, n_polyhedron++)
  //       {
  //         /// loop by all cell vertices
  //         int ncell = *it_polyhedron;
  //         int ivrtx = 0;

  //         vector<int>::iterator it_cell_vrtx;
  //         for ( it_cell_vrtx = vsCellCustom[ncell].vVertices.begin() ; it_cell_vrtx < vsCellCustom[ncell].vVertices.end(); ++it_cell_vrtx, ivrtx++ )
  //         {
  //           if ( *it_cell_vrtx == *it_face_vrtx )
  //             break;
  //         }
  //         vVerticesPair[n_polyhedron] = vsCellCustom[ncell].vVerticesNewnum[ivrtx];
  //       }
  //      vsetGlueVerticies[ vVerticesPair[0] ].insert ( vVerticesPair[1] );
  //      vsetGlueVerticies[ vVerticesPair[1] ].insert ( vVerticesPair[0] );
  //    }
  //   }
  //   /*
  //   if ( vsFaceCustom[ iface ].nMarker > 0 && vsFaceCustom[ iface ].nNeighbors == 2 )
  //   {
  //     vector<int>::iterator it_face_vrtx;

  //     /// loop by all face vertices
  //     for ( it_face_vrtx = vsFaceCustom[ iface ].vVertices.begin();  it_face_vrtx != vsFaceCustom[ iface ].vVertices.end(); ++it_face_vrtx )
  //     {
  //       if ( exlude_some_verticies[*it_face_vrtx] == true)
  //       {
  //         /// loop by neighbor cells
  //         int n_polyhedron = 0;
  //         set<int>::iterator it_polyhedron;
  //         for ( it_polyhedron = vsetPolygonPolyhedron[iface].begin(); it_polyhedron != vsetPolygonPolyhedron[iface].end(); ++it_polyhedron, n_polyhedron++ )
  //         {
  //           /// loop by all cell vertices
  //           int ncell = *it_polyhedron;
  //           int ivrtx = 0;

  //           vector<int>::iterator it_cell_vrtx;
  //           for ( it_cell_vrtx = vsCellCustom[ncell].vVertices.begin() ; it_cell_vrtx < vsCellCustom[ncell].vVertices.end(); ++it_cell_vrtx, ivrtx++ )
  //           {
  //             if ( *it_cell_vrtx == *it_face_vrtx )
  //               break;
  //           }
  //           vVerticesPair[n_polyhedron] = vsCellCustom[ncell].vVerticesNewnum[ivrtx];
  //         }
  //         vsetGlueVerticies[ vVerticesPair[0] ].insert ( vVerticesPair[1] );
  //         vsetGlueVerticies[ vVerticesPair[1] ].insert ( vVerticesPair[0] );
  //       }
  //    }
  //   } */
  // }

  // cout << endl << "Distinguish authenic vertices (might be very slow)" << endl;
  // int n_possible_verticies = vsetGlueVerticies.size();
  // vector<int> v_buf_storage;
  // int total_size_old = 0;
  // int total_size_new = 1;
  // int cycle = 0;
  // while ( total_size_old != total_size_new )
  // {
  //   total_size_old = total_size_new;
  //   total_size_new = 0;
  //   cout << endl << "cycle   :" << cycle << "\t \t Hopefully < 10"; cycle++;

  //   for ( int ivrtx = vsCellCustom[0].nVertices; ivrtx < n_possible_verticies; ivrtx++ )
  //   {
  //     int job_percent = int ( ( 100. * ivrtx ) / ( n_possible_verticies ) );
  //     cout << "\r    " << job_percent << "%";
  //     v_buf_storage.clear();
  //     set<int>::iterator it_set;

  //     for ( it_set = vsetGlueVerticies[ivrtx].begin(); it_set != vsetGlueVerticies[ivrtx].end(); ++it_set )
  //     {
  //       set<int>::iterator it_set_down;

  //       for ( it_set_down = vsetGlueVerticies[ *it_set ].begin(); it_set_down != vsetGlueVerticies[ *it_set ].end(); ++it_set_down )
  //       {
  //         v_buf_storage.push_back ( *it_set_down );
  //       }
  //     }

  //     vector<int>::iterator it_vec;

  //     for ( it_vec = v_buf_storage.begin(); it_vec != v_buf_storage.end(); it_vec++ )
  //     {
  //       vsetGlueVerticies[ivrtx].insert ( *it_vec );
  //     }
  //     total_size_new += vsetGlueVerticies[ivrtx].size();
  //   }
  // }

  // cout << endl << "Renumber vector of stick vertices" << endl;
  // vector<int> vRenumVerticies;
  // vRenumVerticies.resize(n_possible_verticies);

  // /// take first cell
  // for(int ivrtx = 0; ivrtx < vsCellCustom[0].nVertices; ivrtx++)
  // {
  //   vRenumVerticies[ivrtx] = ivrtx;
  // }

  // counter_ = vsCellCustom[0].nVertices;
  // for(int ivrtx = vsCellCustom[0].nVertices; ivrtx < n_possible_verticies; ivrtx++)
  // {
  //   if( *vsetGlueVerticies[ivrtx].begin() == ivrtx )
  //   {
  //     // this vertex is not stick
  //     vRenumVerticies[ivrtx] = counter_;
  //     counter_++;
  //   }
  //   else
  //   {
  //     // this vertex is stick
  //     // set is sorted by c++ defaults, and we take minimum value
  //     vRenumVerticies[ivrtx] = vRenumVerticies[ *vsetGlueVerticies[ivrtx].begin() ];
  //   }
  // }

  // for(int icell = 0; icell < nCells; icell++)
  // {
  //   for(int ivrtx = 0; ivrtx < vsCellCustom[icell].nVertices; ivrtx++)
  //   {
  //     vsCellCustom[icell].vVerticesNewnum[ ivrtx ] = vRenumVerticies[ vsCellCustom[icell].vVerticesNewnum[ ivrtx ] ];
  //   }
  // }
  // int nTotalNodes = counter_;

  // cout << "\t check renumbering consistency " << endl;
  // vector<double> vStatus;
  // vStatus.resize(counter_, false);
  // for ( int icell = 0; icell < nCells; icell++ )
  // {
  //   for ( int ivrtx = 0; ivrtx < vsCellCustom[icell].nVertices; ivrtx++ )
  //     vStatus[ vsCellCustom[icell].vVerticesNewnum[ivrtx] ] = true;
  // }

  // cout << "\t change face numbering" << endl;
  // for(int iface = 0; iface < nFaces; iface++)
  // {
  //   vsFaceCustom[iface].vVerticesNewnum.resize( vsFaceCustom[iface].nVertices, -1);
  //   // we take always [0] - support cell
  //   int icell = vsFaceCustom[iface].vNeighbors[0];
  //   for(int ivrtx = 0; ivrtx < vsFaceCustom[iface].nVertices; ivrtx++)
  //   {
  //     for(int inode = 0; inode < vsCellCustom[icell].nVertices; inode++)
  //     {
  //       if( vsFaceCustom[iface].vVertices[ivrtx] == vsCellCustom[icell].vVertices[inode] )
  //       {
  //         vsFaceCustom[iface].vVerticesNewnum[ivrtx] = vsCellCustom[icell].vVerticesNewnum[inode];
  //         break;
  //       }
  //     }
  //   }
  // }

  // cout << "\t create atoms list" << endl;
  // set<int>::iterator itintset;
  // vector<int> vTempAtoms;
  // vTempAtoms.resize( nTotalNodes, -999 );
  // for (itintset = setIdenticalInternalMarker.begin(); itintset != setIdenticalInternalMarker.end(); ++itintset)
  // {
  //   for(int iface = 0; iface < nFaces; iface++)
  //   {
  //     int ncell = vsFaceCustom[iface].vNeighbors[1];
  //     if( vsFaceCustom[iface].nMarker == *itintset)
  //     {
  //       for(int ivrtx = 0; ivrtx < vsFaceCustom[iface].nVertices; ivrtx++)
  //       {
  //         for(int inode = 0; inode < vsCellCustom[ncell].nVertices; inode++)
  //         {
  //           if( vsFaceCustom[iface].vVertices[ivrtx] == vsCellCustom[ncell].vVertices[inode] )
  //           {
  //             vTempAtoms[ vsCellCustom[ncell].vVerticesNewnum[inode] ] = vsFaceCustom[iface].vVerticesNewnum[ivrtx];
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  // cout << endl << "Change coordanates vector" << endl;
  // vector<vector<double> > vvNewCoordinates;
  // vvNewCoordinates.resize( nTotalNodes, vector<double>(3,0) );
  // for(int icell = 0; icell < nCells; icell++)
  // {
  //   vector<int>::iterator it_old, it_new;
  //   it_new = vsCellCustom[icell].vVerticesNewnum.begin();
  //   for(it_old = vsCellCustom[icell].vVertices.begin(); it_old != vsCellCustom[icell].vVertices.end(); ++it_old, ++it_new)
  //   {
  //     vvNewCoordinates[ *it_new ] = vvVrtxCoords[ *it_old ];
  //   }
  // }

  // vvVrtxCoords.resize(nTotalNodes, vector<double>(3,0));
  // for(int ivrtx = 0; ivrtx < nTotalNodes; ivrtx++)
  // {
  //   vvVrtxCoords[ivrtx] = vvNewCoordinates[ivrtx];
  // }

  // cout << endl << "Unify previous and splitted data" << endl;
  // cout << "Verticies : " << nNodes << "\t \t After splitting : " << nTotalNodes << endl;
  // nNodes = nTotalNodes;
  // for(int icell = 0; icell < nCells; icell++)
  // {
  //   vsCellCustom[icell].vVertices = vsCellCustom[icell].vVerticesNewnum;
  // }

  // for(int iface = 0; iface < nFaces; iface++)
  // {
  //   vsFaceCustom[iface].vVertices = vsFaceCustom[iface].vVerticesNewnum;
  // }

  // vvAtoms.resize(nNodes, vector<int>(2,0) );
  // nAtoms = 0;
  // for(int iatom = 0; iatom < nNodes; iatom++)
  // {
  //   if(vTempAtoms[iatom] >= 0)
  //   {
  //     vvAtoms[nAtoms][0] = iatom;
  //     vvAtoms[nAtoms][1] = vTempAtoms[iatom];
  //     nAtoms++;
  //   }
  // }

  // return;
  // cout << endl << "Smart verticies renumbering" << endl;
  // vector<int> vNodesID;
  // vector<int> vIA;
  // vector<int> vJA;
  // vector<int> vRCM;

  // vIA.push_back(0);
  // for(int ic = 0; ic < nCells; ++ic)
  // {
  //   int n = vsCellCustom[ic].vVertices.size();
  //   for(int iv = 0; iv < n; ++iv)
  //   {
  //     vJA.push_back( vsCellCustom[ic].vVertices[iv] );
  //   }
  //   n += vIA[ic];
  //   vIA.push_back(n);
  // }

  // vRCM.assign(nNodes,-1);
  // pRenum->convert(nCells, nNodes, vIA, vJA, vRCM);

  // for(int icell = 0; icell < nCells; icell++)
  // {
  //   for(int iv = 0; iv < vsCellCustom[icell].vVertices.size(); iv++)
  //     vsCellCustom[icell].vVertices[iv] = vRCM[vsCellCustom[icell].vVertices[iv]];
  // }

  // for(int iface = 0; iface < nFaces; iface++)
  // {
  //   for(int iv = 0; iv < vsFaceCustom[iface].vVertices.size(); iv++)
  //     vsFaceCustom[iface].vVertices[iv] = vRCM[vsFaceCustom[iface].vVertices[iv]];
  // }

  // for(int iatom = 0; iatom < nNodes; iatom++)
  // {
  //   if(vTempAtoms[iatom] >= 0)
  //   {
  //     vvAtoms[nAtoms][0] = vRCM[vvAtoms[nAtoms][0]];
  //     vvAtoms[nAtoms][1] = vRCM[vvAtoms[nAtoms][1]];
  //   }
  // }

  // for(int ivrtx = 0; ivrtx < nTotalNodes; ivrtx++)
  //   vvNewCoordinates[vRCM[ivrtx]] = vvVrtxCoords[ivrtx];

  // for(int ivrtx = 0; ivrtx < nTotalNodes; ivrtx++)
  //   vvVrtxCoords[ivrtx] = vvNewCoordinates[ivrtx];
}


void SimData::handleConnections()
{
  // we always start to count fractures
  int counter = 0;
  for ( int iface = 0; iface < nFaces; iface++ )
  {
    if ( vsFaceCustom[iface].nMarker > 0 && vsFaceCustom[iface].nMarker < 1111110 )
    {
      vsFaceCustom[iface].fluidElement = counter; counter++;
    }

    if ( vsFaceCustom[iface].nMarker < 0 )
      vsFaceCustom[iface].fluidElement = -1;
  }

  // then we count cells
  for(int icell = 0; icell < nCells; icell++)
  {
    vsCellCustom[icell].fluidElement = counter; counter++;
  }

}

#include <random>
void SimData::definePhysicalFacets()
{
  std::size_t n_facets = 0;
  nNeumannFaces = 0;
  nDirichletFaces = 0;
  nDirichletNodes = 0;

  vsPhysicalFacet.resize(nFaces);
  for (int iface = 0; iface < nFaces; iface++)
  {
    if(vsFaceCustom[iface].nMarker < 0)  // probably external face
      for (const auto & conf : config.bcs)
      {
        // std::cout << "key = "<< conf.label << std::endl;

        if( vsFaceCustom[iface].nMarker == conf.label)
        {
          vsPhysicalFacet[n_facets].nface = n_facets;
          vsPhysicalFacet[n_facets].ntype = conf.type;
          vsPhysicalFacet[n_facets].nmark = conf.label;
          vsPhysicalFacet[n_facets].condition = conf.value;
          // std::cout << "doing something boundary "
          //           << conf.label
          //           << std::endl;

          n_facets++;
          if (conf.type == 1)
            nDirichletFaces++;
          else
            nNeumannFaces++;
        }
      }
  }

  nPhysicalFacets = n_facets;


  // int nbnd = vsPhysicalBoundary.size();

  // int nfacets = 0;
  // int nfluid = 0;
  // vsPhysicalFacet.resize(nFaces);
  // for(int iface = 0; iface < nFaces; iface++)
  // {
  //   if( vsFaceCustom[iface].nMarker < 0)
  //   {
  //     for(int i = 0; i < nbnd; i++)
  //     {
  //       if( vsFaceCustom[iface].nMarker == vsPhysicalBoundary[i].nmark)
  //       {
  //         vsPhysicalFacet[nfacets].nface = nfacets;
  //         vsPhysicalFacet[nfacets].ntype = vsPhysicalBoundary[i].ntype;
  //         vsPhysicalFacet[nfacets].nmark = vsPhysicalBoundary[i].nmark;
  //         vsPhysicalFacet[nfacets].vCondition.resize( vsPhysicalBoundary[i].vCondition.size() );
  //         vsPhysicalFacet[nfacets].vCondition = vsPhysicalBoundary[i].vCondition;

  //         vsPhysicalFacet[nfacets].nfluid = -1;
  //         if(vsPhysicalBoundary[i].ntype == 1)
  //         {
  //           nDirichletFaces++;
  //         }
	//   else
  //         {
  //           nNeumannFaces++;
	//     vsPhysicalFacet[nfacets].vCondition[0] = vsPhysicalBoundary[i].vCondition[0];
  //         }
  //         nfacets++;
  //       }
  //     }
  //   }

  //   if ( vsFaceCustom[iface].nMarker > 0 && vsFaceCustom[iface].nMarker < 1111110 )
  //   {
  //     vsPhysicalFacet[nfacets].nface = nfacets;
  //     vsPhysicalFacet[nfacets].ntype = 0;
  //     vsPhysicalFacet[nfacets].nmark = vsFaceCustom[iface].nMarker;
  //     vsPhysicalFacet[nfacets].nfluid = nfluid;
  //     nfacets++;
  //     nfluid++;
  //   }
  // }

  //std::default_random_engine generator;
  //std::default_random_engine generator2;
  //std::normal_distribution<double> distribution(100.0, 20.0);
  //std::normal_distribution<double> distribution2(10.0, 2.0); //mD
  // for ( int iface = 0; iface < nFaces; iface++ )
  // {
  //   vsFaceCustom[iface].aperture = 2.0e-3; //m
  //   vsFaceCustom[iface].conductivity = 0.24e-3 * 0.24e-3 * 0.24e-3 / 12. / 1e-15; //mD.m
  //   if ( vsFaceCustom[iface].nMarker == 1 )
  //   {
  //     vsFaceCustom[iface].aperture = 2.0e-3; //m
  //     vsFaceCustom[iface].conductivity = 0.24e-3 * 0.24e-3 * 0.24e-3 / 12. / 1e-15; //mD.m
  //   }
  // }
}

void SimData::defineStressAndDispOnBoundary()
{
  // PhysicalFace pfFace;
  // int inputBC = vvsBCIn.size();

  // vvsBCOut.resize( inputBC );
  // for(int ibnd = 0; ibnd < inputBC; ibnd++)
  // {
  //   int nfacets = -1;
  //   for ( int iface = 0; iface < nFaces; iface++ )
  //   {
  //     if ( vsFaceCustom[iface].nMarker != 0 ) nfacets++;

  //     if ( vsFaceCustom[iface].nMarker < 0 )
  //     {
  //       for(int i = 0; i < vvsBCIn[ibnd].size(); i++ )
  //       {
  //         if( vsFaceCustom[iface].nMarker == vvsBCIn[ibnd][i].nmark )
  //         {
  //           pfFace.nface = nfacets;
  //           pfFace.ntype = vvsBCIn[ibnd][i].ntype;
  //           pfFace.nmark = vvsBCIn[ibnd][i].nmark;
  //           pfFace.nfluid = -1;
  //           pfFace.axle = vvsBCIn[ibnd][i].axle;

  //           pfFace.condition.clear();
  //           pfFace.vCondition.push_back( vvsBCIn[ibnd][i].vCondition[0] );
  //           // pfFace.condition[i] =;
  //           vvsBCOut[ibnd].push_back( pfFace );
  //         }
  //       }
  //     }
  //   }
  // }
}

void SimData::createSimpleWells()
{
  vector<double> center(3,0);
  /// choise support cell for internal facets
  // FAULT/FRACTURE PART
  for ( int iwell = 0; iwell < nWells; iwell++ )
  {
    int icell = 0;
    for ( int iface = 0; iface < nFaces; iface++ )
    {
      center[0] = vsFaceCustom[iface].center[0];
      center[1] = vsFaceCustom[iface].center[1];
      center[2] = vsFaceCustom[iface].center[2];

      if ( vsFaceCustom[iface].nMarker > 0 )
      {
	if ( abs ( center[0] - vsWell[iwell].vWellCoordinate[0] ) < vsWell[iwell].radius_poisk && abs ( center[1] - vsWell[iwell].vWellCoordinate[1] ) < vsWell[iwell].radius_poisk )
        {
          if(center[0] < -250 || center[0] > 250)
	  {
	    vsWell[iwell].vRadiusPoisk.push_back ( abs ( center[0] - vsWell[iwell].vWellCoordinate[0] ) );
	    vsWell[iwell].vID.push_back ( icell );
	    vsWell[iwell].vWi.push_back ( 100.0 );
	  }
        }
        ++icell;
      }
    }
  }
  // MATRIX PART
  /*
  int n = 0;
  for(int iface = 0; iface < nFaces; iface++)
  {
    if(vsFaceCustom[iface].nMarker > 0) n++;
  }
  for ( int iwell = 0; iwell < nWells; iwell++ )
  {
    int icell = 0;
    for ( int ic = 0; ic < nCells; ic++ )
    {
      center[0] = vsCellCustom[ic].center[0];
      center[1] = vsCellCustom[ic].center[1];
      center[2] = vsCellCustom[ic].center[2];
      if ( abs ( center[0] - vsWell[iwell].vWellCoordinate[0] ) < vsWell[iwell].radius_poisk
	 && (center[2] > vsWell[iwell].vWellCoordinate[2]) && (center[2] < vsWell[iwell].vWellCoordinate[3]) )
      {
        vsWell[iwell].vRadiusPoisk.push_back ( abs ( center[0] - vsWell[iwell].vWellCoordinate[0] ) )  ; //&& (center[2] > vsWell[iwell].vWellCoordinate[2]) && (center[2] < vsWell[iwell].vWellCoordinate[3])
        vsWell[iwell].vID.push_back ( ic + n);
        vsWell[iwell].vWi.push_back ( 10.0 ); //well index [default: 10.0]
      }
    }
    if(vsWell[iwell].vWi.size() < 1)
    {
      cout << "Well definition is wrong" << endl;
    }
  }*/

  for ( int iwell = 0; iwell < nWells; iwell++ )
  {
    vsWell[iwell].datum = 1e16;
    for(int i = 0; i < vsWell[iwell].vID.size(); ++i)
    {
      int icell = vsWell[iwell].vID.size();
      if( fabs(vsCellCustom[icell].center[2]) < vsWell[iwell].datum )
	      vsWell[iwell].datum = abs(vsCellCustom[icell].center[2]);
    }
  }

#if 0
  // limit number of perforations
  for ( int iwell = 0; iwell < nWells; iwell++ )
  {
    vector<double> vrad_;
    vrad_.assign(vsWell[iwell].vRadiusPoisk.size(),0);
    vrad_ = vsWell[iwell].vRadiusPoisk;
    sort(vrad_.begin(), vrad_.end());

    int n_ = std::min(2, int(vrad_.size()));
    for(int j = 0; j < n_; ++j)
    {
      int m_ = 0;
      for(int k = j; k < vsWell[iwell].vRadiusPoisk.size(); ++k)
      {
       cout << vrad_[j] << "\t" << vsWell[iwell].vRadiusPoisk[k] << endl;
       if(abs(vrad_[j] - vsWell[iwell].vRadiusPoisk[k]) < 1e-8 )
       {
	  m_ = k;
          break;
       }
      }
      double rad_ = vsWell[iwell].vRadiusPoisk[j];
      vsWell[iwell].vRadiusPoisk[j] = vsWell[iwell].vRadiusPoisk[m_];
      vsWell[iwell].vRadiusPoisk[m_] = rad_;
      int l_ = vsWell[iwell].vID[j];
      vsWell[iwell].vID[j] = vsWell[iwell].vID[m_];
      vsWell[iwell].vID[m_] = l_;
      vsWell[iwell].vWi[j] = vsWell[iwell].vWi[m_];
    }
    vsWell[iwell].vRadiusPoisk.resize(n_);
    vsWell[iwell].vID.resize(n_);
    vsWell[iwell].vWi.resize(n_);

    for(int k_ = 0; k_ < vsWell[iwell].vRadiusPoisk.size(); ++k_)
    {
      cout << "WELL(" << iwell << "), Perf(" << k_ << ") = "  <<  vsWell[iwell].vID[k_] << endl;
    }
  }
#endif
}
