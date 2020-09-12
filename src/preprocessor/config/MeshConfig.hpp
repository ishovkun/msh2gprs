#pragma once
#include "mesh/CartesianMeshParameters.hpp"

enum class MeshType
{
  cartesian, radial, file
};

struct MeshConfig
{
  MeshType type;
  std::string file;
  mesh::CartesianMeshParameters cartesian;
};
