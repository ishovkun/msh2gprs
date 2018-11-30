#pragma once

#include <vector>

namespace mesh
{

struct Face
{
  std::vector<std::size_t> neighbors;
  int marker = 0;  // internal face by default
  std::size_t index;
  int vtk_id;
};

}
