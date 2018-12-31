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
#include <unordered_set>

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
#include <Well.hpp>


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
  std::vector<double> v_props;
  // poro; perm, perm_x, perm_y, perm_z thc, thc_x, thc_y, thc_z;
  // temp; heat_capacity;

  // biot_plas; biot_flow; biot; young; poisson; density; poron;
  // temperature; pressure; volmult;

  // ref_pres; ref_temp; ref_stres; ref_strain;

  // cohesion; friction; dilation; thermal_expansion; pore_thermal_expansion;
  // vector<double> zmf;
  // vector<double> stress;
};


struct PhysicalFace
{
  int ntype;
  std::size_t nface;
  int nmark;
  int ifracture;
  std::size_t nfluid;
  bool coupled;
  angem::Point<3,double> condition;
  std::vector<std::size_t> neighbor_cells;
  double aperture;
  double conductivity;
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
  mesh::SurfaceMesh<double>           mesh;
  double                              aperture;     // m
  double                              conductivity;  // md-m
};




class SimData
{
public:
  // SimData(const string & inputstream, const SimdataConfig & config);
  SimData(mesh::Mesh & grid, const SimdataConfig & config);
  ~SimData();

  void defineRockProperties();
  void defineEmbeddedFractureProperties();
  void computeCellClipping();
  // void mergeSmallFracCells();
  void definePhysicalFacets();
  void defineStressAndDispOnBoundary();

  void setupWells();
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
  void computeTransBetweenDifferentEfracs();


  void createSimpleWells();

  void extractInternalFaces();
  std::size_t n_default_vars() const;

  // get property from cell->v_props by key
  double get_property(const std::size_t cell,
                      const std::string & key) const;

  angem::Point<3,double> get_permeability(const std::size_t cell) const;
  double get_volume_factor(const std::size_t cell) const;
  void meshFractures();

  bool is_fracture (const int marker)
  {
    const auto it = fracture_face_markers.find(marker);
    if (it != fracture_face_markers.end())
      return true;
    else return false;
  }

  bool is_boundary (const int marker)
  {
    const auto it = boundary_face_markers.find(marker);
    if (it != boundary_face_markers.end())
      return true;
    else return false;
  }

protected:
  void compute_frac_frac_intersection_transes(const std::vector<angem::Point<3,double>>   & verts,
                                              const std::vector<std::vector<std::size_t>> & polys,
                                              const std::vector<int>                      & markers,
                                              FlowData                                    & flow_data) const;
  std::size_t get_flow_element_index(const std::size_t ifrac,
                                     const std::size_t ielement) const;
  void setupSimpleWell(Well & well);
  void setupComplexWell(Well & well);

  // renum * pRenum;

public:
  // user-defined program config defined in json or yaml files
  SimdataConfig config;
  // class that handles mesh (cells, faces, neighbors, face-splitting)
  mesh::Mesh & grid;
  // class that stores dfm grid for vtk output
  mesh::SurfaceMesh<double> dfm_master_grid;

  // container for cell properties (user-defined)
  std::vector<RockProps> vsCellRockProps;
  // ass user-defined names for rock properties
  std::vector<std::string> rockPropNames;

  // stores embedded fracture mechanics data
  vector<EmbeddedFracture> vEfrac;
  // stores flow cell volumes, trances, etc.
  FlowData flow_data;
  // might be different from flow_data if the user requests to mesh
  // embedded fractures independently
  FlowData new_flow_data;

  // stores faces with mechanics neumann and dirichlet boundary conditions
  std::unordered_map<std::size_t, PhysicalFace> boundary_faces;
  // stores faces that represent dfm fractures
  std::unordered_map<std::size_t, PhysicalFace> dfm_faces;
  // number of flow dfm faces (before split)
  std::size_t n_flow_dfm_faces;

  // number of faces with dirichlet mechanics conditions
  std::size_t n_dirichlet_faces;
  // number of faces with neumann mechanics conditions
  std::size_t n_neumann_faces;

  // coupling mechanics and flow
  // these are master DFM faces
  // each dfm frac has 2 sides, but Timur thought It's a good idea to pass just one
  std::vector<std::vector<std::size_t>> gm_cell_to_flow_cell;

  // set of markers for dfm faces (used in is_fracture)
  std::set<int> fracture_face_markers;
  // set of markers for boundary faces (used in is_boundary)
  std::unordered_set<int> boundary_face_markers;

  // old timur wells: rewrute
  // vector<SimpleWell> vsWell;

protected:
  StandardElements * pStdElement;
};
