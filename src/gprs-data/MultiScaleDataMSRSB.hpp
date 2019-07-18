#pragma once

#include "mesh/Mesh.hpp"
#include "LayerDataMSRSB.hpp"
#include "UnionFindWrapper.hpp"
#include "tuple_hash.hpp"
#include <algorithm>  // std::max_element
#include <vector>

namespace multiscale
{

using std::size_t;
using std::vector;

class MultiScaleDataMSRSB
{
 public:
  /* Constructor.
   * takes n_blocks for only a single layer,
   * since multi-level multiscale is a long way
   * down the road. */
  MultiScaleDataMSRSB(mesh::Mesh  & grid,
                      const size_t  n_blocks);
  // get reference to the active layer
  inline
  LayerDataMSRSB & active_layer(){return layers[active_layer_index];}
  const LayerDataMSRSB & active_layer() const {return layers[active_layer_index];}

 protected:

 private:
  // call to metis to obtain partitioning
  void build_partitioning();
  // main method that identifies regions where shape functions exist
  void build_support_regions();
  // find geometric centers of coarse blocks
  void find_centroids();
  // build connection map that stores faces between blocks and their centers
  void build_block_connections();
  // identify ghost blocks, find block face centroids,
  // find block edge centroids
  void build_block_face_data();
  // find ghost coarse blocks
  // Note: ghost cells are like top, bottom, etc.
  // their indices start from layer.n_blocks and go up
  void build_ghost_block_faces();
  // check whether two faces share an edge
  bool share_edge(const mesh::const_face_iterator &face1,
                  const mesh::const_face_iterator &face2);
  // build a structure that diistinguishes between boundary face groups
  // this is done to find ghost blocks
  algorithms::UnionFindWrapper<size_t> build_external_face_disjoint();

  bool is_ghost_block(const size_t block) const;

  void find_block_face_centroids(algorithms::UnionFindWrapper<size_t> & face_disjoint,
                                 std::unordered_map<size_t, size_t>   & map_block_group);

  /* collect vertices from faces that are on block-block interfaces */
  std::unordered_map<std::size_t, std::vector<std::size_t>>
  build_map_vertex_blocks(algorithms::UnionFindWrapper<size_t> & face_disjoint,
                          std::unordered_map<size_t, size_t>   & map_block_group);

  //this method inverts a map obtained in the previous method
  // build a map (triplet of block in ascending order) -> (vertex)
  // this map essentially stores block edges and corners
  std::unordered_map<std::tuple<std::size_t,std::size_t,std::size_t>, std::vector<std::size_t>>
  build_block_edges(const std::unordered_map<std::size_t, std::vector<std::size_t>> & map_block_vertices);

  // find centers of block edges and modify layer.block_faces structure
  void find_block_edge_centroids(
      const
      std::unordered_map<std::tuple<std::size_t,std::size_t,std::size_t>, std::vector<std::size_t>>
      &map_block_vertices);



  // members
  const mesh::Mesh & grid;
  vector<LayerDataMSRSB> layers;
  size_t active_layer_index;
};


// void MultiScaleDataMSRSB::()
// {

// }


}
