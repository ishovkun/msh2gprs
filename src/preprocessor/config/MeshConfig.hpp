#pragma once
#include "mesh/CartesianMeshParameters.hpp"
#include "angem/Rotation.hpp"


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
  // how much padding to add between well nodes and outer domain boundary
  // this value should be > 0.
  double padding_fraction = 1.f;
  int cell_label = 0;
};

struct MeshConfig
{
  MeshType type;
  std::string file;
  mesh::CartesianMeshParameters cartesian;
  RefinementParameters refinement;
  INSIMMeshConfig insim;
  std::vector<angem::Rotation<double>> rotations;
};
