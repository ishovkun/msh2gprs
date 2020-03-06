#pragma once

#include "PreprocessorConfig.hpp"
#include "SimData.hpp"
#include "discretization/flow/DoFNumbering.hpp"
#include <set>

namespace gprs_data {

class DiscreteFractureManager
{
 public:
  DiscreteFractureManager(const std::vector<DiscreteFractureConfig> & config,
                          SimData & data);
  // count the number of dfm faces
  size_t count_dfm_faces() const;
  /* Assign flow properties to dfm control volumes */
  void distribute_properties();
  /* returns true if the marker corresponds to a fracture marker */
  bool is_fracture(const int face_marker) const;
  /* split grid dfm faces for geomechanics */
  void split_faces(mesh::Mesh & grid);
  // combine two vectors
  static std::vector<DiscreteFractureConfig>
  combine_configs(const std::vector<DiscreteFractureConfig> & config1,
                  const std::vector<DiscreteFractureConfig> & config2);
  // return vector of split fracture face markers
  std::vector<int> get_face_markers() const;
  // build dfm surface grid for vtk output
  mesh::SurfaceMesh<double> build_dfm_grid(const mesh::Mesh & grid) const;
  // map dfm surface grid to flow dofs
  std::vector<size_t> map_dfm_grid_to_flow_dofs(const mesh::Mesh & grid,
                                                const discretization::DoFNumbering & dofs) const;

 protected:
  void build_dfm_markers_set_();

  const std::vector<DiscreteFractureConfig> & m_config;
  mesh::Mesh & m_grid;
  SimData & m_data;
  std::set<int> m_dfm_markers;  // set of dfm markers
};

}  // end namespace gprs_data
