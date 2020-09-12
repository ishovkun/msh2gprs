#pragma once
#include "angem/Point.hpp"
#include <vector>

namespace mesh {

struct CartesianMeshParameters
{
  std::vector<double> dx = {1.0};
  std::vector<double> dy = {1.0};
  std::vector<double> dz = {1.0};
  angem::Point<3, double> origin;
};

}  // end namespace mesh
