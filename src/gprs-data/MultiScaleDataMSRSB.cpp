#include "MultiScaleDataMSRSB.hpp"
#include "MetisInterface.hpp"


namespace multiscale {

using Point = angem::Point<3,double>;

MultiScaleDataMSRSB::MultiScaleDataMSRSB(mesh::Mesh  & grid,
                                         const size_t  n_blocks)
    :
    grid(grid),
    active_layer_index(0)
{
  auto & layer = layers.emplace_back();
  layer.index = 0;
  layer.n_blocks = n_blocks;
  layer.n_cells = grid.n_cells();

  build_partitioning();
  build_support_regions();
}


void MultiScaleDataMSRSB::build_partitioning()
{
    std::cout << "building connection map" << std::endl;
    auto & layer = active_layer();

    for (auto it = grid.begin_faces(); it != grid.end_faces(); ++it)
    {
      const auto & neighbors = it.neighbors();
      if (neighbors.size() == 2)  // not a boundary face
        layer.cell_connections.insert_connection( neighbors[0], neighbors[1] );
    }

    layer.partitioning = multiscale::MetisInterface<hash_algorithms::empty>
        ::build_partitioning(layer.cell_connections, layer.n_blocks, layer.n_cells);
}


void MultiScaleDataMSRSB::build_support_regions()
{
  find_centroids();
  auto block_face_centroids = find_block_face_centroids();

}


void MultiScaleDataMSRSB::find_centroids()
{
  auto & layer = active_layer();
  layer.block_centroids.resize(layer.n_blocks);
  vector<size_t> n_cells_per_block(layer.n_blocks);

  for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
  {
    const size_t block = layer.partitioning[cell.index()];
    layer.block_centroids[block] += cell.center();
    n_cells_per_block[block]++;
  }

  for (std::size_t block=0; block<layer.n_blocks; ++block)
    layer.block_centroids[block] /= n_cells_per_block[block];
}


void MultiScaleDataMSRSB::find_block_face_centroids()
{
  auto & layer = active_layer();
  vector<Point> face_centers(grid.n_faces());
  vector<size_t> n_sub_faces(grid.n_faces(), 0);

  for (const auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    const auto & neighbors = face.neighbors();
    if (neighbors.size() == 2)  // not a boundary face
    {
      const std::size_t i1 = layer.partitioning[neighbors[0]];
      const std::size_t i2 = layer.partitioning[neighbors[1]];;
      if (i1 != i2)
      {
        const std::size_t connection_index = layer.p_block_connections->connection_index(i1, i2);
        face_centers[connection_index] += face.center;
        n_sub_faces[connection_index]++;
      }

    }
  }
  for (std::size_t i=0; i<layer.p_cell_connections->size(); ++i)
    face_centers[i] /= n_sub_faces[i];
  return face_centers;

}

}  // end namespace
