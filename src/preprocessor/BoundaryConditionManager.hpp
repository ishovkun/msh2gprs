#pragma once
#include "PreprocessorConfig.hpp"
#include "SimData.hpp"
#include "muparser/muParser.h" // parser for user-defined expressions for reservoir data

namespace gprs_data {

/* Manager for mechanical boundary conditions */
class BoundaryConditionManager
{
 public:
  BoundaryConditionManager(const std::vector<BCConfig> & face_config,
                           const std::vector<BCConfig> & node_config,
                           SimData & data);

 private:
  // build boundary conditions and save them into simdata
  void build_boundary_conditions_();
  // create a config for a neumann face
  void process_neumann_face_(const mesh::Face & face, const size_t config_index);
  // process a dirichlet face (put node configs in tmp storage)
  void process_dirichlet_face_(const mesh::Face & face, const size_t config_index);
  // take each non-empty entry in _node_to_config, create a config for it and save into simdata for output
  void create_dirichlet_data_();
  // parse location muparser expressions and find boundary labels
  void find_faces_from_expressions_();
  //
  std::vector<mu::Parser> create_location_parsers_(const std::vector<size_t> & configs);
  void create_value_parsers_(const std::vector<BCConfig> & config,
                             std::vector<std::array<mu::Parser,3>> & parsers);
  // --------------- Variables --------------//
  const std::vector<BCConfig> _face_config;
  const std::vector<BCConfig> _node_config;
  SimData & _data;

  // map node to face configs
  // if no config, it's not a boudnary node
  // if more than one config - it's a node on edge shared by several dirichlet boundaries
  std::vector<std::vector<size_t>> _node_to_config;
  std::array<double,4> _variables;  // variables for muparser (X,Y,Z,nan)
  std::vector<std::array<mu::Parser,3>> _face_value_parsers;
  std::vector<std::array<mu::Parser,3>> _node_value_parsers;
};

}  // end namespace gprs_data
