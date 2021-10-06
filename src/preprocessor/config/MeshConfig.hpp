#pragma once
#include "mesh/CartesianMeshParameters.hpp"

enum class MeshType
{
  cartesian, radial, file, insim
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

struct INSIMMeshConfig
{
  double minimum_thickness = 1.f;  // handy if we want pseudo 2d mesh
};

struct MeshConfig
{
  MeshType type;
  std::string file;
  mesh::CartesianMeshParameters cartesian;
  RefinementParameters refinement;
  INSIMMeshConfig insim;
};
