#pragma once

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <map>

#include "element.hpp"
#include "renum.hpp"
#include "transes.hpp"
// #include "GElement.hpp"

#include "angem/Point.hpp"
#include "angem/PolyGroup.hpp"
#include "angem/Collisions.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "mesh/Mesh.hpp"
#include "SimdataConfig.hpp"
#include <Well.hpp>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>
#include <set>
#include <unordered_set>


namespace gprs_data
{

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


// struct SimpleWell
// {
//   vector<double> vRadiusPoisk;
//   vector<double> vWellCoordinate;
//   vector<int> vID;
//   vector<int> vWi;
//   double datum;
//   string Type;
//   double radius_poisk;
// };


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




/* This class is the body of the preprocessor.
 * It literally contains all the capabilities of the software
 * described in the wiki
 */
class SimData
{
public:
  // SimData(const string & inputstream, const SimdataConfig & config);
  SimData(mesh::Mesh & grid, const SimdataConfig & config);
  // destructor
  ~SimData();

  // distribute user-defined properties over cells
  void defineRockProperties();
  // determine which cells are touched by embedded fractures
  // and set mechanical SDA properties
  void defineEmbeddedFractureProperties();
  // determine geometry of intersection of embedded fractures with the mesh
  void computeCellClipping();
  // merge small section elemenents of edfm fractures with larger neighbors (flow only)
  // the implementation is for the old architecture, needs a rewrite
  void mergeSmallFracCells();
  // create cartesian mesh within edfm fractures for flow computation
  void meshFractures();
  // find boundary and dfm faces in the original mesh
  void definePhysicalFacets();
  // define the cells occupied by well and compute well productivities
  void setupWells();
  // split dfm faces for geomechanics
  void splitInternalFaces();
  // determine which flow volumes correspond to which mechanics cells
  void handleConnections();
  // compute flow data without edfm fracs (reservoir and dfm only)
  void computeReservoirTransmissibilities();
  // high-level method that combines everything related to embedded frac flow
  void handleEmbeddedFractures();
  // compute flow data withich a single edfm fracture
  void computeFracFracTran(const std::size_t                 frac,
                           const EmbeddedFracture          & efrac,
                           const mesh::SurfaceMesh<double> & mesh,
                           flow::FlowData                  & frac_flow_data);
  // compute flow data for edfm fractures
  void computeEDFMTransmissibilities(const std::vector<angem::PolyGroup<double>> & splits,
                                     const int   frac_ind);
  // compute flow data between two edfm fractures --may be old impl
  void computeTransEfracIntersection();
  // number of default variables (such as cell x,y,z) for rock properties
  std::size_t n_default_vars() const;
  // get property from cell->v_props by key
  double get_property(const std::size_t cell,
                      const std::string & key) const;
  // wrapper around get_property that aborts if no perm data available
  angem::Point<3,double> get_permeability(const std::size_t cell) const;
  // wrapper around get_property that aborts if no perm data available
  double get_volume_factor(const std::size_t cell) const;

  // helper: check if face is a fracture
  bool is_fracture (const int marker)
  {
    const auto it = fracture_face_markers.find(marker);
    if (it != fracture_face_markers.end())
      return true;
    else return false;
  }

  // helper: check if face is a boundary face
  bool is_boundary (const int marker)
  {
    const auto it = boundary_face_markers.find(marker);
    if (it != boundary_face_markers.end())
      return true;
    else return false;
  }

protected:
  // compute flow data between two edfm fracs
  void compute_frac_frac_intersection_transes(const std::vector<angem::Point<3,double>>   & verts,
                                              const std::vector<std::vector<std::size_t>> & polys,
                                              const std::vector<int>                      & markers,
                                              flow::FlowData                              & flow_data) const;
  // get flow volume index of an edfm element
  // std::size_t get_flow_element_index(const std::size_t ifrac,
  //                                    const std::size_t ielement) const;
  // create a well that occupies a single cell in z direction
  void setupSimpleWell(Well & well);
  // create a complex well that occupies multiple cells and is arbitrarily-oriented
  void setupComplexWell(Well & well);
  // compute productivities of all well segments
  void computeWellIndex(Well & well);
  // get dimensions of a cell bounding box
  angem::Point<3,double> get_dx_dy_dz(const std::size_t icell) const;

  // is given flow element an embedded fracture
  bool is_embedded_fracture(const std::size_t flow_element_index) const;
  // is given flow element a discrete fracture
  bool is_discrete_fracture(const std::size_t flow_element_index) const;
  // is given flow element a reservoir cell
  bool is_reservoir_element(const std::size_t flow_element_index) const;
  // return global flow index of an ielement element of embedded fracture i
  std::size_t efrac_flow_index(const std::size_t ifrac,
                               const std::size_t ielement) const;
  std::size_t res_cell_flow_index(const std::size_t icell) const {return n_flow_dfm_faces + icell;}

  void apply_projection_edfm(const std::size_t ifrac,     // embedded frac index
                             const std::size_t ielement,  // frac element index
                             const std::size_t icell);    // reservoir cell index

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
  flow::FlowData flow_data;
  // might be different from flow_data if the user requests to mesh
  // embedded fractures independently
  flow::FlowData new_flow_data;
  // wells
  std::vector<Well> wells;

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
  // i'm not sure if it's even used
  StandardElements * pStdElement;
  // class that performs vertex renumbering  after dfm split for faster computation
  // renum * pRenum;
};

}
