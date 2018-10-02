#pragma once

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <map>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>
#include <set>
using namespace std;

#include "element.hpp"
#include "renum.hpp"

#include "Point.hpp"
#include "PolyGroup.hpp"
#include "Collisions.hpp"
#include <SimdataConfig.hpp>

struct GmConstraint
{
  double dPenalty;
  int nVertices;
  vector<int> vIndexes;
  vector<int> vVertices;
  vector<vector<double> > vvCoefficient;
};


struct RockProps
{
  int model;
  std::vector<double> v_props;
  double poro;
  double perm, perm_x, perm_y, perm_z;
  double thc, thc_x, thc_y, thc_z;
  double temp;
  double heat_capacity;

  double biot_plas;
  double biot_flow;
  double biot;
  double young;
  double poisson;
  double density;
  double poron;

  double temperature;
  double pressure;
  double volmult;

  double ref_pres;
  double ref_temp;
  double ref_stres;
  double ref_strain;

  double cohesion;
  double friction;
  double dilation;
  double thermal_expansion;
  double pore_thermal_expansion;

  vector<double> zmf;
  vector<double> stress;
};

struct Gelement
{
  int nMarker;
  int nVertices;
  int nNeighbors;

  int vtkIndex, formIndex;
  int fluidElement;

  vector<std::size_t> vVertices;
  vector<std::size_t> vVerticesSorted;
  vector<int> vVerticesNewnum;

  vector<std::size_t> vNeighbors;

  double thickness, center_distance;
  // vector<double> vCenter;
  // vector<double> vNormal;
  angem::Point<3,double> center;
  angem::Point<3,double> normal;
  double aperture, conductivity;
};

struct PhysicalFace
{
  int ntype;
  int nface;
  int nmark;
  int axle;
  int nfluid;
  // vector<double> vCondition;
  angem::Point<3,double> condition;
};


struct SimpleWell
{
  vector<double> vRadiusPoisk;
  vector<double> vWellCoordinate;
  vector<int> vID;
  vector<int> vWi;
  double datum;
  string Type;
  double radius_poisk;
};


struct EmbeddedFracture
{
  std::vector<std::size_t>            cells;
  std::vector<angem::Point<3,double>> points;
  std::vector<double>                 dip;
  std::vector<double>                 strike;
  double                              cohesion;
  double                              friction_angle;
  double                              dilation_angle;
  // cells -> points
  // these two entries represent mesh within the frac
  std::vector<angem::Point<3,double>>  vVertices;
  // cells -> vertex indiced
  std::vector<std::vector<std::size_t>> vIndices;
  // std::vector<ScratchData>              vSplits;
};




class SimData
{
public:
  SimData(const string & inputstream, const SimdataConfig & config);
  ~SimData();
  // void readSetupValues();
  void readTotalData();
  void readTotalTemp();

  void readGmshFile();
  void convertGmsh2Sim();

  void initilizeBoundaryConditions();

  void defineRockProperties();
  void defineEmbeddedFractureProperties();
  void computeCellClipping();
  void definePhysicalFacets();
  void defineStressAndDispOnBoundary();

  void splitInternalFaces();

  void handleConnections();

  void createSimpleWells();

  void extractInternalFaces();
  std::size_t n_default_vars();

  // get property from cell->v_props by key
  double get_property(const std::size_t cell,
                      const std::string & key) const;

  angem::Point<3,double> get_permeability(const std::size_t cell) const;

protected:
   void methodElementCenter(int nelem, vector<Gelement> &vsElement);
   void methodFaceNormalVector(int nelem, vector<Gelement> &vsElement);
   void methodChangeFacesNormalVector();

   void methodRandomRockProperties();
   double createLognormalDistribution(double E, double S);

   int checkReservedBoundaryName(int nmarker)
   {
     if( nmarker > 1000000 ) return(-1);
     return(1);
   }
   renum * pRenum;

public:
  int corner_cell;
  bool dual_media;
  vector<double> grade_total, temp_total;
  int nNodes;

  double dNotNumber;

  int nBndNodes;
  vector<vector<double> > vvBndFaceNodesCoor;
  vector<vector<int> >    vvBndFaceNodes;
  vector<int>    vBndFaceCode;

  vector<vector<double> > vvInputCoorNodes;
  vector<vector<int> >    vvElementNodes;
  vector<int> vElementCode;

  string instream;
  string outstream;

  // Internal Data
  std::size_t nVertices;
  double maxVrtxCoordsX;
  double maxVrtxCoordsY;
  double maxVrtxCoordsZ;
  double minVrtxCoordsX;
  double minVrtxCoordsY;
  double minVrtxCoordsZ;
  vector< angem::Point<3,double> > vvVrtxCoords;
  vector<vector<double> > vvVrtxDisp;
  vector<int> vConstraintVertex;

  int nCells;
  vector<Gelement> vsCellCustom;
  vector<RockProps> vsCellRockProps;
  vector<EmbeddedFracture> vEfrac;

  int nAtoms;
  vector<vector<int> > vvAtoms;

  int nExternalBoundaryFaces;
  int nInternalBoundaryFaces;
  std::size_t nFaces;
  vector<Gelement> vsFaceCustom;

  set<int> setIdenticalExternalMarker;
  set<int> setIdenticalInternalMarker;

  vector<double> vIdenticalInternalFacetPerm;
  vector<double> vIdenticalInternalFacetAperture;
  vector<double> vIdenticalInternalFacetFFpermMult;
  vector<double> vIdenticalInternalFacetFMpermMult;

  // polygons and polyhedrons
  vector<set<int> > vsetPolyhedronPolygon;
  vector<set<int> > vsetPolygonPolyhedron;

  int nDirichletFaces;
  int nDirichletNodes;
  int nNeumannFaces;

  int nPhysicalFacets;
  vector<PhysicalFace> vsPhysicalFacet;
  vector<PhysicalFace> vsPhysicalBoundary;

  vector<vector<PhysicalFace> > vvsBCIn;
  vector<vector<PhysicalFace> > vvsBCOut;

  //wells
  int nWells;
  vector<SimpleWell> vsWell;

  double Sxx, Syy, Szz, Syz, Sxz, Sxy;

  SimdataConfig config;


protected:
  StandardElements * pStdElement;
  void computePropertyMaps();

  struct tokens: std::ctype<char>
  {
    tokens(): std::ctype<char>(get_table()) {}

    static std::ctype_base::mask const* get_table()
    {
        typedef std::ctype<char> cctype;
        static const cctype::mask *const_rc= cctype::classic_table();

        static cctype::mask rc[cctype::table_size];

        rc[' '] = std::ctype_base::space;
        rc['\t'] = std::ctype_base::space;
        rc['/'] = std::ctype_base::space;
        return &rc[0];
   }
  };

  struct wordtokens: std::ctype<char>
  {
    wordtokens(): std::ctype<char>(get_table()) {}

    static std::ctype_base::mask const* get_table()
    {
        typedef std::ctype<char> cctype;
        static const cctype::mask *const_rc= cctype::classic_table();

        static cctype::mask rc[cctype::table_size];

        rc['_'] = std::ctype_base::space;
        rc['.'] = std::ctype_base::space;
        return &rc[0];
   }
  };

  vector<int> vPointPass;
  vector<double> vPointCoord;
  vector<vector<double> > vvPlate;
};
