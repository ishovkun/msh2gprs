#pragma once

#include "PreprocessorConfig.hpp"
// #include "ControlVolumeData.hpp"
// #include "ConnectionData.hpp"
#include "SimData.hpp"

namespace gprs_data {

class EmbeddedFractureManager
{
 public:
  EmbeddedFractureManager(std::vector<EmbeddedFractureConfig> &config,
                          const EDFMMethod edfm_method,
                          SimData & data);
  void split_cells();
  // generate DiscreteFractureConfig object that form due to
  // cell splitting
  std::vector<DiscreteFractureConfig> generate_dfm_config();
  // true if face marker belongs to an edfm fracture after splitting cells
  bool is_fracture(const int face_marker) const;
  // distribute SDA properties
  void distribute_mechanical_properties();
  // map SDA cells to edfm control volumes
  void map_mechanics_to_control_volumes();
  // return vector of split fracture face markers
  std::vector<int> get_face_markers() const;

 private:
  bool find_edfm_cells_(angem::Polygon<double> & fracture, std::vector<size_t> & cells);
  // split internal grid cells due to intersection with embedded fracture
  void split_cells_(angem::Polygon<double> & fracture, std::vector<size_t> & cells, const int face_marker);
  // find the maximum face marker of the grid
  int find_maximum_face_marker_() const;
  // ------------------ Variables -----------------
  // non-const cause we move fractures to avoid collision with vertices
  std::vector<EmbeddedFractureConfig> &config;
  // simple or pedfm
  EDFMMethod m_method;
  // will be filled
  SimData & data;
  mesh::Mesh & m_grid;
  std::set<int> m_edfm_markers;
  // std::unordered_map<size_t> m_edfm_faces;
};

}  // end namespace gprs_data
