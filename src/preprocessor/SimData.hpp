#pragma once

#include "discretization/ControlVolumeData.hpp"
#include "discretization/ConnectionData.hpp"
#include "mesh/Mesh.hpp"
#include "Well.hpp"
#include "angem/Tensor2.hpp"
#include <unordered_map>

namespace gprs_data {

/* Data of a single dfm fracture grid face */
struct DiscreteFractureFace
{
  int    marker;                        // fracture marker
  size_t cv_index;                      // index of the control volume
  bool   coupled;                       // coupling with geomechanics
  double aperture;                      // hydraulic aperture of the fracture [m]
  double conductivity;                  // hydraulic conductivity of dfm fracture [m·md]
  std::vector<double> custom_flow_data;
};

struct EmbeddedFractureMechanicalProperties
{
  std::vector<std::size_t>  cells;            // cells that the fracture crosses
  std::vector<angem::Point<3,double>> points; // points in the frac plane within the intersected cells
  std::vector<double> dip;                    // fracture dip angle in a cell [°]
  std::vector<double> strike;                 // fracture strike angle in a cell [°]
  double cohesion;                            // fracture cohesive strength [bar]
  double friction_angle;                      // fracture friction angle [°]
  double dilation_angle;                      // fracture dilation angle [°]
  // double aperture;                            // hydfraulic aperture [m]
  // double conductivity;                        // hydraulic conductivity [md-m]
  mesh::SurfaceMesh<double> mesh;             // combined grid discretization of all embedded fractures
};

struct SimData
{
  mesh::Mesh grid;
  // cell properties
  // ----------------------- Reservoir cells ------------------------ //
  std::vector<std::string> property_names;
  std::vector<std::vector<double>> cell_properties;
  std::array<double,9> permeability_keys;  // permeability key indices in cell_properties
  size_t porosity_key_index;               // porosity key index in cell_properties
  std::vector<size_t> output_flow_properties;
  std::vector<size_t> cell_cv_indices;  // flow control indices of reservoir cells
  // ----------------------- DFM ------------------------ //
  std::unordered_map<size_t,DiscreteFractureFace> dfm_faces;
  // grid comprised of dfm faces
  mesh::SurfaceMesh<double> dfm_grid;
  // ---------------------- EDFM ------------------------ //
  std::vector<EmbeddedFractureMechanicalProperties> sda_data;
  // ----------------------- Flow data ---------------------- //
  std::vector<discretization::ControlVolumeData> cv_data;
  std::vector<discretization::ConnectionData> flow_connection_data;
  // ----------------------- Well data ---------------------- //
  std::vector<Well> wells;  // vector of well properties
  angem::PointSet<3,double> well_vertices;  // set of well coordinatees: used for vtk output.
  // vector of well segments: indices of well coordinate points. used for vtk output.
  std::vector<std::pair<std::size_t,std::size_t>> well_vertex_indices;
  // --------------------- Methods --------------------------------- //
  angem::Tensor2<3,double> get_permeability(const std::size_t cell) const
  {
    assert(cell < cell_properties[permeability_keys[0]].size());
    angem::Tensor2<3,double> K;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        K(i, j) = (permeability_keys[3*i+j] >= 0) ?
                  cell_properties[permeability_keys[3*i+j]][cell] : 0;
    return K;
  }

  double get_porosity(const std::size_t cell) const
  {
    assert(cell < property_names.size());
    assert(porosity_key_index > 0 && porosity_key_index < property_names.size());
    return cell_properties[porosity_key_index][cell];
  }

  size_t n_dfm_faces() const { return dfm_faces.size(); }
};


}  // end namespace gprs_data
