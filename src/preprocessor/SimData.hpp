#pragma once

#include "mesh/Mesh.hpp"

namespace gprs_data {

struct SimData
{
  mesh::Mesh grid;
  // cell properties
  std::vector<std::string> property_names;
  std::vector<std::vector<double>> cell_properties;
};


}  // end namespace gprs_data
