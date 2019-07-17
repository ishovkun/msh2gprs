#pragma once

#include "ConnectionMap.hpp"
#include "angem/Point.hpp"

namespace multiscale
{

struct BlockFace
{
  angem::Point<3,double> center;
  std::vector<std::size_t> edging_blocks;
  std::vector<angem::Point<3,double>> edge_centers;
};

struct LayerDataMSRSB
{
  std::size_t index;  // index of the layer
  // number of coarse blocks
  std::size_t n_blocks;
  std::size_t n_cells;
  // block (coarse scale) index in each (fine) cell
  std::vector<std::size_t> partitioning;  // size = n_fine_cells
  // inverse of the partitioning
  //  block (coarse cell) -> list of fine cells in it
  std::vector<std::vector<std::size_t>> cells_in_block;
  // connections on a lower level (fine-scale connections)
  // std::shared_ptr<PureConnectionMap> p_cell_connections;
  PureConnectionMap cell_connections;
  // coarse connections
  // Note: not empty on the last level (filled explicitly during partitioning)
  PureConnectionMap block_connections;
  // faces between blocks and blocks and ghost cells
  hash_algorithms::ConnectionMap<BlockFace> block_faces_data;

  // ghost cell data (external faces)
  std::size_t n_ghost_cells;
  std::vector<std::vector<std::size_t>> adjacent_ghost_cells;
  std::vector<angem::Point<3,double>> ghost_cell_directions;
  std::vector<std::string> ghost_cell_names;

  // std::shared_ptr< std::vector<angem::Point<3,double>> > p_cell_centroids;
  std::vector<angem::Point<3,double>>  block_centroids;
  // ---------------- Flowe     ( fvm ) stuff ----------------------------- //
  std::shared_ptr< std::vector<angem::Point<3,double>> > p_cell_centroids_flow;
  std::shared_ptr< std::vector<angem::Point<3,double>> > p_block_centroids_flow;
  // indices of centroid cells
  std::vector<std::size_t> block_centroid_cell_indices;
  // ---------------- mechanics ( fem ) stuff ----------------------------- //
  // indices of centroid vertices
  std::vector<std::size_t> block_centroid_vertex_indices;
  // support region
  std::vector<std::unordered_set<std::size_t>> support_boundary_cells;
  std::vector<std::unordered_set<std::size_t>> support_internal_cells;

  // i impose dirichet here
  std::vector<std::vector<std::size_t>> support_boundary_nodes;
  // i store shape function values in these nodes
  std::vector<std::vector<std::size_t>>   support_internal_nodes;
  // std::vector<std::vector<double>> shape_function_values;
  // std::vector<std::unordered_map<std::size_t, double> > shape_function_values;
};

}  // end namespace
