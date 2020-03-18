#pragma once
#include "PreprocessorConfig.hpp"
#include "SimData.hpp"

namespace gprs_data {

/* Manager for mechanical boundary conditions */
class BoundaryConditionManager
{
 public:
  BoundaryConditionManager(const std::vector<BCConfig> & face_config,
                           const std::vector<BCNodeConfig> & node_config,
                           SimData & data);

 private:
  // build boundary conditions and save them into simdata
  void build_boundary_conditions_();
  // create a config for a neumann face
  void process_neumann_face_(const BCConfig & conf, const size_t face_index);
  // process a dirichlet face (put node configs in tmp storage)
  void process_dirichlet_face_(const mesh::Face & face, const size_t config_index);
  // take each non-empty entry in _node_to_config, create a config for it and save into simdata for output
  void create_dirichlet_data_();
  // --------------- Variables --------------//
  const std::vector<BCConfig> _face_config;
  const std::vector<BCNodeConfig> _node_config;
  SimData & _data;

  // map node to face configs
  // if no config, it's not a boudnary node
  // if more than one config - it's a node on edge shared by several dirichlet boundaries
  std::vector<std::vector<size_t>> _node_to_config;
};

}  // end namespace gprs_data
