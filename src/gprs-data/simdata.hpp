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

#include "element.hpp"
#include "renum.hpp"
#include "transes.hpp"
#include "GElement.hpp"

#include "angem/Point.hpp"
#include "angem/PolyGroup.hpp"
#include "angem/Collisions.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/Mesh.hpp"
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
  mesh::SurfaceMesh<double>            mesh;
  std::vector<angem::Point<3,double>>  vVertices;
  // cells -> vertex indiced (fracture polygons)
  std::vector<std::vector<std::size_t>> vIndices;
  double                               aperture;     // m
  double                               conductivity;  // md-m
};




class SimData
{
public:
  // SimData(const string & inputstream, const SimdataConfig & config);
  SimData(mesh::Mesh & grid, const SimdataConfig & config);
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
  // void mergeSmallFracCells();
  void definePhysicalFacets();
  void defineStressAndDispOnBoundary();

  void splitInternalFaces();

  void handleConnections();
  void computeReservoirTransmissibilities();
  void computeFracFracTran(const std::size_t                 frac,
                           const EmbeddedFracture          & efrac,
                           const mesh::SurfaceMesh<double> & mesh,
                           FlowData                        & frac_flow_data);
  void computeEDFMTransmissibilities(const std::vector<angem::PolyGroup<double>> & splits,
                                     const int   frac_ind);
  void computeInterEDFMTransmissibilities();
  // void computeTransBetweenDifferentEfracs();


  void createSimpleWells();

  void extractInternalFaces();
  std::size_t n_default_vars() const;

  // get property from cell->v_props by key
  double get_property(const std::size_t cell,
                      const std::string & key) const;

  angem::Point<3,double> get_permeability(const std::size_t cell) const;
  double get_volume_factor(const std::size_t cell) const;
  void meshFractures();

protected:
   void methodElementCenter(int nelem, vector<Gelement> &vsElement);
   void methodFaceNormalVector(int nelem, vector<Gelement> &vsElement);
   void methodChangeFacesNormalVector();

   void methodRandomRockProperties();
   double createLognormalDistribution(double E, double S);

  void handle_edfm_face_intersection(const std::size_t ifrac,
                                     const std::size_t jfrac,
                                     const std::vector<std::size_t> & icells,
                                     const std::vector<std::size_t> & jcells);

  void compute_frac_frac_intersection_transes(const std::vector<angem::Point<3,double>>   & verts,
                                              const std::vector<std::vector<std::size_t>> & polys,
                                              const std::vector<int>                      & markers,
                                              FlowData                                    & flow_data) const;
  // std::size_t get_flow_element_index(const std::size_t ifrac,
  //                                    const std::size_t ielement) const;


   int checkReservedBoundaryName(int nmarker)
   {
     if( nmarker > 1000000 ) return(-1);
     return(1);
   }
   renum * pRenum;

public:
  mesh::Mesh & grid;
  int corner_cell;
  bool dual_media;
  vector<double> grade_total, temp_total;
  // int nNodes;
  std::size_t nNodes;

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

  // int nCells;
  // vector<Gelement> vsCellCustom;
  std::vector<RockProps> vsCellRockProps;
  std::vector<std::string> rockPropNames;

  vector<EmbeddedFracture> vEfrac;
  FlowData flow_data;
  FlowData new_flow_data;

  int nAtoms;
  vector<vector<int> > vvAtoms;

  int nExternalBoundaryFaces;
  int nDFMFracs;
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
        // static const cctype::mask *const_rc= cctype::classic_table();

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
        // static const cctype::mask *const_rc= cctype::classic_table();

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
