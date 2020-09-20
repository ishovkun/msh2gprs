#pragma once
#include "mesh/CartesianMeshParameters.hpp"

enum class MeshType
{
  cartesian, radial, file
};

enum class RefinementType
{
  aspect_ratio, none
};

struct RefinementParameters
{
  size_t max_level = 10;
  double aspect_ratio;
  RefinementType type = RefinementType::none;
};

struct MeshConfig
{
  MeshType type;
  std::string file;
  mesh::CartesianMeshParameters cartesian;
  RefinementParameters refinement;
};
