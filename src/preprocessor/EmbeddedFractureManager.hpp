#pragma once

#include "PreprocessorConfig.hpp"
#include "discretization/DoFNumbering.hpp"
#include "SimData.hpp"

namespace gprs_data {

class EmbeddedFractureManager
{
 public:
  /* Constructor */
  EmbeddedFractureManager(std::vector<EmbeddedFractureConfig> &config,
                          const EDFMMethod edfm_method,
                          const double min_dist_to_node,
                          SimData & data);
  void split_cells();
  // generate DiscreteFractureConfig object that form due to
  // cell splitting
  std::vector<DiscreteFractureConfig> generate_dfm_config();
  // true if face marker belongs to an edfm fracture after splitting cells
  bool is_fracture(const int face_marker) const;
  // distribute SDA properties
  void distribute_mechanical_properties();
  // return vector of split fracture face markers
  std::vector<int> get_face_markers() const;
  // build edfm surface grid for vtk output
  void build_edfm_grid(const discretization::DoFNumbering & dofs);
  // map SDA cells to edfm control volumes
  // do it only after coarsening the grid
  // and distribute mechanical properties
  void map_mechanics_to_control_volumes(const discretization::DoFNumbering & dofs);

 private:
  bool find_edfm_cells_(angem::Polygon<double> & fracture, std::vector<size_t> & cells);
  void find_edfm_cells_and_faces_();
  // split internal grid cells due to intersection with embedded fracture
  void split_cells_(angem::Polygon<double> & fracture, std::vector<size_t> & cells, const int face_marker);
  // find the maximum face marker of the grid
  int find_maximum_face_marker_() const;
  // wrapper around m_marker_config;
  size_t fracture_index_(const int face_marker) const;
  // ------------------ Variables -----------------
  // non-const cause we move fractures to avoid collision with vertices
  std::vector<EmbeddedFractureConfig> &config;
  // simple or pedfm
  EDFMMethod m_method;
  // minimum distance from fracture to vertex relative to the cell size
  const double _min_dist_to_node;
  // will be filled
  SimData & m_data;
  mesh::Mesh & m_grid;
  // marker to config index
  std::map<int,size_t> m_marker_config;
};

}  // end namespace gprs_data
