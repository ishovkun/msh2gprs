#pragma once

#include <vector>

namespace multiscale
{

struct MultiScaleOutputData
{
  // size = n_fine_cells
  std::vector<std::size_t> partitioning;
  // cells that constitute boundaries of each support region
  std::vector<std::vector<std::size_t>> support_boundary_cells;
  // cells that dconstitute the internals of each support region
  // this includes the cells inside the coarse block
  std::vector<std::vector<std::size_t>> support_internal_cells;
};

}
