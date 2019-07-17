#include "MultiScaleDataMSRSB.hpp"
#include "MetisInterface.hpp"


namespace multiscale {

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



}  // end namespace
