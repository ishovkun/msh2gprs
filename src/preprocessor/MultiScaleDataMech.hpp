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
  virtual void fill_output_model(MultiScaleOutputData & model, const int layer_index = 0) const override;

 protected:
  // find coarse block corners
  // void find_block_corners(const algorithms::UnionFindWrapper<size_t> & face_disjoint,
  //                         const std::unordered_map<size_t, size_t>   & map_block_group,
  //                         const std::vector<std::unordered_set<std::size_t>> & cell_block_neighbors);
  void find_block_corners(const std::unordered_map<size_t, size_t> & map_boundary_face_ghost_block,
                          const std::vector<std::unordered_set<std::size_t>> & cell_block_neighbors);
  void find_block_corners2(const std::unordered_map<size_t, size_t> & map_boundary_face_ghost_block,
                          const std::vector<std::unordered_set<std::size_t>> & cell_block_neighbors);

  // build vector cell -> list of blocks it is a neighbors of
  std::vector<std::unordered_set<std::size_t>>
  // build_cell_block_neighbors(const algorithms::UnionFindWrapper<size_t> & face_disjoint,
  //                            const std::unordered_map<size_t, size_t>   & map_block_group) const;
  build_cell_block_neighbors(const std::unordered_map<size_t, size_t> & map_boundary_face_ghost_block) const;
  // build a container for boundary vertices
  void build_boundary_nodes(const std::unordered_map<size_t, size_t> & map_boundary_face_ghost_block);

  // const mesh::Mesh & grid;
  // std::unordered_map<size_t, std::unordered_set<size_t>> map_coarse_node_blocks;
  std::vector<std::vector<size_t>> coarse_node_blocks;
  // std::vector<std::unordered_set<size_t>> coarse_node_blocks;
  // std::unordered_map<size_t,size_t> map_boundary_face_ghost_block;

};

}  // end namespca
