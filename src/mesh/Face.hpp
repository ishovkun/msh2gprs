#pragma once

#include <vector>

namespace mesh
{
static const int DEFAULT_FACE_MARKER = -1;

struct Face
{
  std::vector<std::size_t> neighbor_cells;
  int marker = DEFAULT_FACE_MARKER;  // internal face by default
  std::size_t index;
  std::size_t master_face_index;  // index before split
  int vtk_id;
  std::vector<std::size_t> vertices;
};

}
