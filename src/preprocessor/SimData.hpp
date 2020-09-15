#pragma once

#include "discretization/flow/DoFNumbering.hpp"           // provides discretization::DoFNumbering
#include "discretization/flow/ControlVolumeData.hpp"      // provides discretization::ControlVolumeData
#include "discretization/flow/ConnectionData.hpp"         // provides discretization::ConnectionData
#include "discretization/mechanics/FiniteElementData.hpp" // provides FiniteElementData
#include "intersections/GridIntersectionSearcher.hpp"     // provides GridIntersectionSearcher
#include "mesh/Mesh.hpp"                                  // provides mesh::Mesh
#include "Well.hpp"                                       // provides Well
#include "MultiScaleOutputData.hpp"                       // provides multiscale::MultiScaleOutputData
#include "angem/Tensor2.hpp"                              // provides angem::Tensor2
#include <unordered_map>

namespace gprs_data {

/* Data of a single dfm fracture grid face */
struct DiscreteFractureFace
{
  int    marker;                        // fracture marker
  bool   coupled;                       // coupling with geomechanics
  double aperture;                      // hydraulic aperture of the fracture [m]
  double conductivity;                  // hydraulic conductivity of dfm fracture [m·md]
  std::vector<double> custom_flow_data;
  size_t region;  // fracture region (used by AD-GPRS for property tables)
};

struct EmbeddedFractureMechanicalProperties
{
  std::vector<std::size_t>  cells;            // unsplit mechanics cells that the fracture crosses
  std::vector<std::vector<size_t>> faces;     // split face indices
  std::vector<angem::Point<3,double>> points; // points in the frac plane within the intersected cells
  std::vector<double> dip;                    // fracture dip angle in a cell [°]
  std::vector<double> strike;                 // fracture strike angle in a cell [°]
  double cohesion;                            // fracture cohesive strength [bar]
  double friction_angle;                      // fracture friction angle [°]
  double dilation_angle;                      // fracture dilation angle [°]
  double conductivity;                        // hydraulic conductivity of edfm fracture [m·md]
  mesh::SurfaceMesh<double> mesh;             // combined grid discretization of all embedded fractures
  size_t region;                              // fracture region (used by AD-GPRS for property tables)
};

struct BoundaryConstraintData
{
  std::vector<size_t> nodes;
  std::vector<size_t> components;
  double penalty;
};

struct SimData
{
  mesh::Mesh grid;  // active grid that has all the manipulations on
  mesh::Mesh geomechanics_grid;  // mechanics grid copied before cell edfm splitting
  // cell properties
  // ----------------------- Reservoir cells ------------------------ //
  std::vector<std::string> property_names;
  std::vector<std::vector<double>> cell_properties;
  std::array<int,6> permeability_keys = {-1, -1, -1, -1, -1, -1}; // permeability key indices in cell_properties
  size_t porosity_key_index = std::numeric_limits<size_t>::max(); // porosity key index in cell_properties
  size_t vol_mult_index = std::numeric_limits<size_t>::max();     // volume mult index in cell_properties
  std::vector<size_t> output_flow_properties;                     // indices of flow property keywords
  std::vector<size_t> output_mech_properties;                     // indices of mech property keywords
  // ----------------------- DFM ------------------------ //
  // std::unordered_map<size_t,DiscreteFractureFace> dfm_faces;
  std::map<size_t,DiscreteFractureFace> dfm_faces;
  // grid comprised of dfm faces
  // mesh::SurfaceMesh<double> dfm_flow_grid, dfm_mech_grid;
  mesh::SurfaceMesh<double> fracture_grid;
  std::vector<size_t> dfm_cell_mapping;  // for postprocessor output  vtk_cell -> flow dof
  std::vector<angem::Point<3,double>> grid_vertices_after_face_split;
  std::vector<std::vector<size_t>> grid_cells_after_face_split;
  // map parent vertex -> vector of child vertices after face splitting for mech DFM
  std::unordered_map<size_t, std::vector<size_t>> parent_to_child_vertices;
  // ---------------------- EDFM ------------------------ //
  std::vector<EmbeddedFractureMechanicalProperties> sda_data;
  // std::unordered_map<size_t,size_t> face_to_fracture;
  // mesh:: SurfaceMesh<double> edfm_grid;     // for vtk output
  std::vector<size_t> edfm_cell_mapping;    // for postprocessor output  vtk_cell -> dof
  // std::unordered_set<int> edfm_grid_labels;
  // ----------------------- Flow data ---------------------- //
  std::vector<discretization::ControlVolumeData> cv_data;
  std::vector<discretization::ConnectionData> flow_connection_data;
  // ----------------------- Well data ---------------------- //
  std::vector<Well> wells;  // vector of well properties
  angem::PointSet<3,double> well_vertices;  // set of well coordinatees: used for vtk output.
  // vector of well segments: indices of well coordinate points. used for vtk output.
  std::vector<std::pair<std::size_t,std::size_t>> well_vertex_indices;
  // =========================== GEOMECHANICS ================= //
  bool has_mechanics = false;
  std::shared_ptr<discretization::DoFNumbering> mech_numbering;  // mech cell and face numbering
  std::shared_ptr<discretization::DoFNumbering> flow_numbering;  // flow dof numbering
  std::vector<bool> coupling;                                    // if grid cells are coupled
  // ----------------------- FEM data  ---------------------- //
  using FEMData = discretization::FiniteElementData;
  std::vector<FEMData> fe_cell_data;  // fe values and gradients for grid cells
  std::vector<FEMData> fe_face_data;  // fe values and gradients for grid faces
  std::vector<std::vector<FEMData>> fe_frac_data;  // fe values and gradients of cells in face qpoints
  // ----------------------- Boundary conditions ------------ //
  std::vector<size_t> neumann_face_indices;  // indices of neumann faces
  std::vector<angem::Point<3,double>> neumann_face_traction;  // values of neuman bc's
  std::array<std::vector<size_t>, 3> dirichlet_indices;       // indices of dirichlet vertices
  std::array<std::vector<double>, 3> dirichlet_values;        // dirichlet values in vertices
  std::vector<BoundaryConstraintData> boundary_constraints;   // penalized constrained vertex groups
  // ----------------------- Multiscale ------------ //
  multiscale::MultiScaleOutputData ms_mech_data;
  // ----------------------- Other ---------------------- //
  std::vector<std::vector<size_t>> gmcell_to_flowcells; // Geomechanics cell -> Flow cells in each geomech cell.
  std::vector<std::vector<size_t>> gmcell_to_SDA_flowcells; // Geomechanics cell in EDFM -> Flow cells in each geomech cell in EDFM.
  // Helps to search iintersections fast
  std::unique_ptr<GridIntersectionSearcher> grid_searcher;
  // --------------------- Methods --------------------------------- //
  angem::Tensor2<3,double> get_permeability(const std::size_t cell) const
  {
    assert(cell < cell_properties[permeability_keys[0]].size());
    angem::Tensor2<3,double> K;
    for (int i = 0; i < 3; i++)
        K(i, i) = (permeability_keys[i] >= 0) ?
                  cell_properties[permeability_keys[i]][cell] : 0;
    return K;
  }

  double get_porosity(const std::size_t cell) const
  {
    assert(porosity_key_index < property_names.size());
    assert ( cell < cell_properties[porosity_key_index].size() );
    return cell_properties[porosity_key_index][cell];
  }

  double get_volume_mult(const size_t cell) const
  {
    if (vol_mult_index == std::numeric_limits<size_t>::max())
      return 1;
    else return cell_properties[vol_mult_index][cell];
  }

  size_t n_dfm_faces() const { return dfm_faces.size(); }
};


}  // end namespace gprs_data
