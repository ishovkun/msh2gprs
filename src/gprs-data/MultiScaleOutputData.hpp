#pragma once

#include <unordered_set>
#include <vector>

namespace multiscale
{

struct MultiScaleOutputData
{
  size_t n_blocks;
  //  size = n_fine_cells
  std::vector<std::size_t> partitioning;
  std::vector<std::size_t> centroids;
  // cells that constitute boundaries of each support region
  std::vector<std::unordered_set<std::size_t>> support_boundary_cells;
  // cells that dconstitute the internals of each support region
  // this includes the cells inside the coarse block
  std::vector<std::unordered_set<std::size_t>> support_internal_cells;
};

}
