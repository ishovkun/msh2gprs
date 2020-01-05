#pragma once

#include "PreprocessorConfig.hpp"
// #include "ControlVolumeData.hpp"
// #include "ConnectionData.hpp"
#include "SimData.hpp"

namespace gprs_data {

class EmbeddedFractureManager
{
 public:
  EmbeddedFractureManager(const std::vector<EmbeddedFractureConfig> &config,
                          const EDFMMethod edfm_method,
                          SimData & data);
  void split_cells();
  // generate DiscreteFractureConfig object that form due to
  // cell splitting
  std::vector<DiscreteFractureConfig> generate_dfm_config();
  // true if face marker belongs to an edfm fracture after splitting cells
  bool is_fracture(const int face_marker) const;

  // extract cv data pertaining to edfm fractures from the mixed assembly
  std::vector<discretization::ControlVolumeData>
  extract_control_volume_data(const std::vector<discretization::ControlVolumeData> & mixed_cv_data,
                              const size_t n_dfm_faces);
  // compute edfm transmissibility from cut-cell transmissibility
  void extract_flow_data(const std::vector<discretization::ControlVolumeData> & mixed_cv_data,
                         const std::vector<discretization::ConnectionData> & mixed_connection_data,
                         const size_t n_dfm_faces, const size_t n_cells);

 private:
  bool find_edfm_cells_(angem::Polygon<double> & fracture,
                        std::vector<size_t> & cells) const;
  // split internal grid cells due to intersection with
  // embedded fracture
  void split_cells_(angem::Polygon<double> & fracture,
                    std::vector<size_t> & cells,
                    const int face_marker);
  // find the maximum face marker of the grid
  int find_maximum_face_marker_() const;
  // ------------------ Variables -----------------
  const std::vector<EmbeddedFractureConfig> &config;
  // simple or pedfm
  EDFMMethod m_method;
  // will be filled
  SimData & data;
  mesh::Mesh & m_grid;
  std::set<int> m_edfm_markers;
  // std::unordered_set<size_t> m_edfm_faces;
};

}  // end namespace gprs_data
