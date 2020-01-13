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
  void create_properties();

 private:
  void create_neumann_faces_(const std::unordered_map<size_t,size_t> & face_to_frac);
  std::unordered_map<size_t,size_t> create_boundary_faces_() const;
  // --------------- Variables --------------//
  const std::vector<BCConfig> m_face_config;
  const std::vector<BCNodeConfig> m_node_config;
  SimData & m_data;
};

}  // end namespace gprs_data
