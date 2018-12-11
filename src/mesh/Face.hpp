#pragma once

#include <vector>

namespace mesh
{

struct Face
{
  std::vector<std::size_t> neighbors;
  int marker = 0;  // internal face by default
  std::size_t index;
  std::size_t old_index;  // index before split
  int vtk_id;
  std::vector<std::size_t> ordered_indices;
};

}
