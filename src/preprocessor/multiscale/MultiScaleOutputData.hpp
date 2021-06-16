#pragma once

#include <unordered_set>
#include <vector>

namespace multiscale
{

struct MultiScaleOutputData
{
  // if true, expert into cells, else export vertex data
  bool cell_data = true;
  size_t n_coarse;
  //  size = n_fine_cells
  std::vector<std::size_t> partitioning;  // it's always cell data cause I said so
  std::vector<std::size_t> centroids;
  // cells that constitute boundaries of each support region
  std::vector<std::vector<std::size_t>> support_boundary;
  // cells that dconstitute the internals of each support region
  // this includes the cells inside the coarse block
  std::vector<std::vector<std::size_t>> support_internal;
};

}
