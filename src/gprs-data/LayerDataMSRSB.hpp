#pragma once

#include "ConnectionMap.hpp"
#include "mesh/SurfaceMesh.hpp"
#include "angem/Point.hpp"

#include <unordered_set>

namespace multiscale
{

struct BlockFace
{
  angem::Point<3,double> center;
  std::vector<std::size_t> edging_blocks;
  std::vector<angem::Point<3,double>> edge_centers;
  std::size_t n_cell_faces = 0;
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
  // coarse block centers
  std::vector<angem::Point<3,double>>  block_centroids;
  // face connections between coarse elements
  PureConnectionMap block_internal_connections;
  // Comprises:
  // (1) internal connections between blocks
  // (2) blocks and ghost cells
  // (3) edge block-block and block-ghost connections
  // Note: ghost cells are like top, bottom, etc.
  // their indices start from layer.n_blocks and go up
  hash_algorithms::ConnectionMap<BlockFace> block_faces;

  // ghost cell data (external faces)
  // std::size_t n_ghost_cells;
  // std::vector<std::vector<std::size_t>> adjacent_ghost_cells;
  // std::vector<angem::Point<3,double>> ghost_cell_directions;
  // std::vector<std::string> ghost_cell_names;
  // ---------------- Flow      ( fvm ) stuff ----------------------------- //
  // std::shared_ptr< std::vector<angem::Point<3,double>> > p_cell_centroids_flow;
  // std::shared_ptr< std::vector<angem::Point<3,double>> > p_block_centroids_flow;
  // indices of centroid cells
  // std::vector<std::size_t> block_centroid_cell_indices;
  // ---------------- Mechanics ( fem ) stuff ----------------------------- //
  // indices of centroid vertices
  std::vector<std::size_t> coarse_node_vertices;
  // support region
  std::vector<std::unordered_set<std::size_t>> support_boundary_cells;
  std::vector<std::unordered_set<std::size_t>> support_internal_cells;

  // i impose dirichet here
  // std::vector<std::vector<std::size_t>> support_boundary_nodes;
  // i store shape function values in these nodes
  // std::vector<std::vector<std::size_t>>   support_internal_nodes;
  // std::vector<std::vector<double>> shape_function_values;
  // std::vector<std::unordered_map<std::size_t, double> > shape_function_values;

  // for vtk output
  // mesh::SurfaceMesh<double> support_bounding_surface;
};

}  // end namespace
