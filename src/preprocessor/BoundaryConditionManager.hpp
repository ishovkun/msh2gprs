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
  const std::vector<BCConfig> m_face_config;
  const std::vector<BCNodeConfig> m_node_config;
  SimData & m_data;
};

}  // end namespace gprs_data
