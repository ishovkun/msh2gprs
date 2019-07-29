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
  // coarse block geometric centers
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

  // support region data
  std::vector<std::unordered_set<std::size_t>> support_boundary;
  std::vector<std::unordered_set<std::size_t>> support_internal;
  // centroid to a vertex or cell
  std::vector<std::size_t> coarse_to_fine;
};

}  // end namespace
