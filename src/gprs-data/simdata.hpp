#pragma once

#include <cstdlib>
#include <string>
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
#include "MultiScaleOutputData.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>
#include <set>
#include <unordered_set>


namespace gprs_data
{

// TODO: replace with just a vector of properties
struct RockProps
{
  std::vector<double> v_props;
};


/* structure hold data about dfm fractures and
 * faces with mechanical boundary conditions */
struct PhysicalFace
{
  //  dirichlet or neumann
  int ntype;
  // face index
  std::size_t nface;
  // gmsh face marker
  int nmark;
  // index of dfm fracture the frace belongs to
  int ifracture;
  // index of fluid control volume
  std::size_t nfluid;
  // whether to map fracture face to reservoir mechanical cell
  bool coupled;
  // boundary condition components
  angem::Point<3,double> condition;
  // faces that the face belongs to (1 or 2)
  std::vector<std::size_t> neighbor_cells;
  // hydraulic aperture of dfm fracture
  double aperture;
  // hydraulic conductivity of dfm fracture
  double conductivity;
};


/* contrains data for embedded fracture geometry,
 * mechanical and flow parameters */
struct EmbeddedFracture
{
  // cells that the fracture crosses
  std::vector<std::size_t>  cells;
  // points in the frac plane within the intersected cells
  std::vector<angem::Point<3,double>> points;
  // for mechanics: fracture dip angle in a cell
  std::vector<double> dip;
  // for mechanics: fracture strike angle in a cell
  std::vector<double> strike;
  // mechanical paramters
  double                    cohesion;   // fracture cohesive strength
  double                    friction_angle; // fracture friction angle
  double                    dilation_angle; // fracture dilation angle
  // combined grid discretization of all embedded fractures
  mesh::SurfaceMesh<double> mesh;
  double                    aperture;   // hydfraulic aperture [m]
  double                    conductivity; // hydraulic conductivity [md-m]
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

  //  distribute user-defined properties over cells
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
  //  high-level method that combines everything related to embedded frac flow
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

  // Multiscale
  void build_multiscale_data();

protected:
  // number of default variables (such as cell x,y,z) for rock properties
  std::size_t n_default_vars() const;
  // get property from cell->v_props by key
  double get_property(const std::size_t cell,
                      const std::string & key) const;
  // wrapper around get_property that aborts if no perm data available
  angem::Point<3,double> get_permeability(const std::size_t cell) const;
  // wrapper around get_property that aborts if no perm data available
  double get_volume_factor(const std::size_t cell) const;
  // compute flow data between two edfm fracs
  void compute_frac_frac_intersection_transes(const std::vector<angem::Point<3,double>>   & verts,
                                              const std::vector<std::vector<std::size_t>> & polys,
                                              const std::vector<int>                      & markers,
                                              flow::FlowData                              & flow_data) const;
  // get flow volume index of an edfm element
  // create a well that occupies a single cell in z direction
  void setupSimpleWell(Well & well);
  // create a complex well that occupies multiple cells and is arbitrarily-oriented
  void setupComplexWell(Well & well);
  // compute productivities of all well segments
  void computeWellIndex(Well & well);
  // get dimensions of a cell bounding box
  angem::Point<3,double> get_dx_dy_dz(const std::size_t icell) const;

  // is given flow element an embedded fracture
  // Params [in]
  // flow_element_index: index of the flow control volume
  //
  // Note:
  // edfm control volumes span from 0 to n_edfm_intersections
  // reservoir control volumes space from n_edfm to n_edfm + n_cells
  // dfm control volumes span from n_edfm + n_cells to n_edfm + n_cells + n_dfm
  bool is_embedded_fracture(const std::size_t flow_element_index) const;
  // whether given flow element a discrete fracture
  // Params [in]
  // flow_element_index: index of the flow control volume
  //
  // Note:
  // edfm control volumes span from 0 to n_edfm_intersections
  // reservoir control volumes space from n_edfm to n_edfm + n_cells
  // dfm control volumes span from n_edfm + n_cells to n_edfm + n_cells + n_dfm
  bool is_discrete_fracture(const std::size_t flow_element_index) const;
  // is given flow element a reservoir cell
  //
  // Note:
  // edfm control volumes span from 0 to n_edfm_intersections
  // reservoir control volumes space from n_edfm to n_edfm + n_cells
  // dfm control volumes span from n_edfm + n_cells to n_edfm + n_cells + n_dfm
  bool is_reservoir_element(const std::size_t flow_element_index) const;
  // return global flow index of an ielement element of embedded fracture i
  std::size_t efrac_flow_index(const std::size_t ifrac,
                               const std::size_t ielement) const;
  std::size_t res_cell_flow_index(const std::size_t icell) const
  {return n_flow_dfm_faces + icell;}

  // connect embedded fractures to cells in a physical way
  // Params [in]
  // ifrac: embedded frac index ndex of embedded fracture
  // ielement: frac element index  / fracture element index
  // icell: reservoir cell index
  // split: result of dissecting a cell with a fracture polygon
  void apply_projection_edfm(const std::size_t                ifrac,
                             const std::size_t                ielement,
                             const std::size_t                icell,
                             const angem::PolyGroup<double> & split);
  // face selection method for pedfm
  // pedfm projects fracture onto the smaller faces of the
  // intersected cell
  // Params [in]
  // cell: iterator pointing to the intersecting cell
  // split: result of dissecting a cell with a fracture polygon
  std::vector<mesh::face_const_iterator>
  pedfm_select_faces(const mesh::cell_iterator      & cell,
                     const angem::PolyGroup<double> & split) const;

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
  // pureply for visulalization // vtk_segments;
  angem::PointSet<3,double> well_vertices;
  std::vector<std::pair<std::size_t,std::size_t>> well_vertex_indices;

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

  // multiscale
  multiscale::MultiScaleOutputData ms_flow_data;
  multiscale::MultiScaleOutputData ms_mech_data;

  // different from partitioning cause of fracturess and wells
  //  std::vector<std::size_t> fluid_partitioning;

protected:
  // it might be used in older timur's version for 2nd order elements but not
  // in this version
  StandardElements * pStdElement;
  // class that performs vertex renumbering  after dfm split for
  // linear solver operation
  // renum * pRenum;
};

}
