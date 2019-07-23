#pragma once

#include "MultiScaleDataMSRSB.hpp"
// #include "LayerDataMSRSB.hpp"
#include "MultiScaleOutputData.hpp"
#include "mesh/Mesh.hpp"
#include <algorithm>  // std::max_element
#include <vector>

namespace multiscale {

using std::size_t;
using std::vector;

class MultiScaleDataMech : public MultiScaleDataMSRSB
{
 public:
  MultiScaleDataMech(mesh::Mesh  & grid, const size_t  n_blocks);
  virtual void build_data() override;
  // virtual void fill_output_model(MultiScaleOutputData & model, const int layer_index = 0) const;

 protected:
  // find coarse block corners
  void find_block_corners(const algorithms::UnionFindWrapper<size_t> & face_disjoint,
                          const std::unordered_map<size_t, size_t>   & map_block_group,
                          const std::vector<std::unordered_set<std::size_t>> & cell_block_neighbors);

  // build vector cell -> list of blocks it is a neighbors of
  std::vector<std::unordered_set<std::size_t>>
  build_cell_block_neighbors(const algorithms::UnionFindWrapper<size_t> & face_disjoint,
                             const std::unordered_map<size_t, size_t>   & map_block_group) const;
  // const mesh::Mesh & grid;
  // size_t active_layer_index;

};

}  // end namespca
