#pragma once

#include "mesh/Mesh.hpp"
#include "LayerDataMSRSB.hpp"
#include "MultiScaleOutputData.hpp"
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
  MultiScaleDataMSRSB(mesh::Mesh  & grid, const size_t  n_blocks);
  // main method. that's when the fun happens
  virtual void build_data();
  virtual void fill_output_model(MultiScaleOutputData & model, const int layer_index = 0) const;

 protected:
  // get reference to the active layer
  inline
  LayerDataMSRSB & active_layer(){return layers[active_layer_index];}
  const LayerDataMSRSB & active_layer() const {return layers[active_layer_index];}

  // call to metis to obtain partitioning
  void build_partitioning();
  // build inverted partitioning block -> cells
  void build_cells_in_block();
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
  bool share_edge(const mesh::Face &face1,
                  const mesh::Face &face2) const;
  // build a structure that diistinguishes between boundary face groups
  // this is done to find ghost blocks
  algorithms::UnionFindWrapper<size_t> build_external_face_disjoint() const;
  // wrapper around a previous function
  std::unordered_map<size_t,size_t> build_map_face_to_ghost_cell() const;


  bool is_ghost_block(const size_t block) const;

  void find_block_face_centroids(const std::unordered_map<size_t,size_t> & map_boundary_face_ghost_block);

  /* collect vertices from faces that are on block-block interfaces */
  std::unordered_map<std::size_t, std::vector<std::size_t>>
  build_map_vertex_blocks(const std::unordered_map<size_t,size_t> & map_boundary_face_ghost_block);

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

  // build support region for a block
  void build_support_region(const std::size_t block);

  // mark cells that lay on the boundary of the support region
  // (intersect with the bounding shape)
  // Input:
  // block: coarse element for which we are building the support region
  // neighbor: coarse block in which the support boundary is located
  // bounding_shape: a triangle that cuts through the cells in the neighbor block
  // those find cells are marked as boundary cells
  void build_support_region_boundary(const std::size_t block,
                                     const std::size_t neighbor,
                                     const angem::Shape<double> & bounding_shape);
  // find a cell that's definitely outside the support region for the current block
  angem::Point<3,double> find_point_outside_support_region(const std::size_t block);

  // build the internal points of the block support region
  void build_support_internal_cells(const std::size_t block,
                                    const mesh::SurfaceMesh<double>& bounding_surface);

  //  attributes
  const mesh::Mesh & grid;
  vector<LayerDataMSRSB> layers;
  size_t active_layer_index;

 private:
  mutable std::unordered_map<size_t, std::string> debug_ghost_cell_names;
  void debug_make_ghost_cell_names(const algorithms::UnionFindWrapper<size_t> & face_disjoint,
                                   std::unordered_map<size_t, size_t>   & map_block_group) const;
};


// void MultiScaleDataMSRSB::()
// {

// }


}
