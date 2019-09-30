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
  MultiScaleDataMech(mesh::Mesh  & grid,
                     const std::array<size_t,3> &  n_blocks,
                     const PartitioningMethod method,
                     const size_t elimination_level);
  virtual void build_data() override;
  virtual void fill_output_model(MultiScaleOutputData & model, const int layer_index = 0) const override;

 protected:
  // find coarse block corners
  // void find_block_corners(const algorithms::UnionFindWrapper<size_t> & face_disjoint,
  //                         const std::unordered_map<size_t, size_t>   & map_block_group,
  //                         const std::vector<std::unordered_set<std::size_t>> & cell_block_neighbors);
  void merge_close_blocks();
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
  // new method for building the boundary nodes after vertex repartitioning
  void find_boundary_nodes(const std::unordered_map<size_t, size_t> & map_boundary_face_ghost_block);

  // find cell neighbours of coarse vertices
  // return vector[coarse_node]
  // map cell -> level of neighboring
  // cells that contain the vertex are level-0 neighbors
  // neighbors of those cells are level-2 neighbors
  std::vector<std::unordered_map<size_t,size_t>> find_node_neighboring_cells(const size_t level = 0) const;

  // 0 means no elimination
  // 1 eliminate coarse node if in the same cell as another coarse node
  const size_t elimination_level;
  // std::unordered_map<size_t, std::unordered_set<size_t>> map_coarse_node_blocks;
  std::vector<std::vector<size_t>> coarse_node_blocks;
};

}  // end namespca
