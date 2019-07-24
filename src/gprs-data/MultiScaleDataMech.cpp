#include "MultiScaleDataMech.hpp"
#include <unordered_set>
#include <unordered_map>
#include <map>

namespace multiscale
{

using std::unordered_map;
using std::unordered_set;
using std::map;

MultiScaleDataMech::MultiScaleDataMech(mesh::Mesh  & grid, const size_t  n_blocks)
    :
    MultiScaleDataMSRSB(grid, n_blocks)
{}


void MultiScaleDataMech::build_data()
{
  build_partitioning();

  algorithms::UnionFindWrapper<size_t> face_disjoint = build_external_face_disjoint();
  size_t ghost_block = active_layer().n_blocks;
  std::unordered_map<size_t, size_t> map_block_group;
  for ( const auto &it_face: face_disjoint.items() )
    if (map_block_group.find(it_face.first) == map_block_group.end())
      map_block_group.insert({ it_face.first, ghost_block++});

  const auto cell_block_neighbors = build_cell_block_neighbors(face_disjoint, map_block_group);
  // const auto map_vertex_blocks = build_map_vertex_blocks(face_disjoint, map_block_group);
  // const auto map_block_edge_vertices = build_block_edges(map_vertex_blocks);

  find_block_corners(face_disjoint, map_block_group, cell_block_neighbors);
  build_boundary_nodes();
}


std::vector<std::unordered_set<std::size_t>>
MultiScaleDataMech::build_cell_block_neighbors(const algorithms::UnionFindWrapper<size_t> & face_disjoint,
                                               const std::unordered_map<size_t, size_t>   & map_block_group) const
{
  const auto & layer = active_layer();
  vector<unordered_set<size_t>> cell_block_neighbors(grid.n_cells());

  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    const auto & neighbors = face.neighbors();
    const size_t cell1 = neighbors[0];
    const size_t block1 = layer.partitioning[cell1];
    if (neighbors.size() == 2)
    {
      const size_t cell2 = neighbors[1];
      const size_t block2 = layer.partitioning[cell2];
      if (block1 != block2)
      {
        cell_block_neighbors[cell1].insert(block1);
        cell_block_neighbors[cell1].insert(block2);
        cell_block_neighbors[cell2].insert(block1);
        cell_block_neighbors[cell2].insert(block2);
      }
    }
    else // if (neighbors.size() == 1)
    {
      const size_t face_group = face_disjoint.group(face.index());
      const size_t block2 = map_block_group.find( face_group )->second;
      cell_block_neighbors[cell1].insert(block1);
      cell_block_neighbors[cell1].insert(block2);
    }
  }
  return cell_block_neighbors;
}


void MultiScaleDataMech::
find_block_corners(const algorithms::UnionFindWrapper<size_t>         & face_disjoint,
                   const std::unordered_map<size_t, size_t>           & map_block_group,
                   const std::vector<std::unordered_set<std::size_t>> & cell_block_neighbors)
{
  auto & layer = active_layer();
  // given a list of cells that have more than two neighboring blocks
  // find corners of the blocks.
  // To do this we loop over the faces of the cell, figure out which faces touch
  // other blocks, and then find a common vertex among those vertices
  // unordered_set<size_t> global_corners;
   for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
     if (cell_block_neighbors[cell.index()].size() > 2)  // contains block corner
     {
       // don't want a hash_map here
       map<size_t, unordered_set<size_t>> vertex_blocks;
       map<size_t, unordered_set<size_t>> vertex_faces;

       // a corner is a node all of whose adjacent faces touch different blocks
       for (const auto & face : cell.faces())
       {
         const auto & neighbors = face.neighbors();
         const size_t cell1 = neighbors[0];
         const size_t block1 = layer.partitioning[cell1];
         if (neighbors.size() == 2)
         {
           const size_t block2 = layer.partitioning[neighbors[1]];
           if (block1 != block2)
             for (const auto & vertex : face.vertex_indices())
             {
               vertex_blocks[vertex].insert(block1);
               vertex_blocks[vertex].insert(block2);
               vertex_faces[vertex].insert(face.index());
             }
         }
         else // if (neighbors.size() == 1)
         {
           const size_t face_group = face_disjoint.group(face.index());
           const size_t block2 = map_block_group.find( face_group )->second;
           for (const auto & vertex : face.vertex_indices())
           {
             vertex_blocks[vertex].insert(block1);
             vertex_blocks[vertex].insert(block2);
             vertex_faces[vertex].insert(face.index());
           }
         }
       }

       for (auto & it : vertex_blocks)
         if (it.second.size() > 2)
           if (vertex_faces[it.first].size() > 2)
           { // it's a corner!
             auto glob_it = map_coarse_node_blocks.find(it.first);
             if (map_coarse_node_blocks.find(it.first) == map_coarse_node_blocks.end())
             {
               std::unordered_set<size_t> blocks;
               for (const size_t b : it.second) blocks.insert(b);
               map_coarse_node_blocks.insert({ it.first, std::move(blocks)});
             }
             else
             {
               auto & blocks = glob_it->second;
               for (const size_t b : it.second) blocks.insert(b);
             }
           }
     }

   // save coarse nodes
   layer.coarse_to_fine.resize(map_coarse_node_blocks.size());
   size_t index = 0;
   for (const auto & it : map_coarse_node_blocks)
   {
     layer.coarse_to_fine[index] = it.first;
     index++;
   }
}

void MultiScaleDataMech::fill_output_model(MultiScaleOutputData & model, const int layer_index) const
{
  const auto & layer = layers[layer_index];

  model.cell_data = false;
  model.n_coarse = layer.coarse_to_fine.size();

  // partitioning
  model.partitioning.resize(layer.partitioning.size());
  std::copy(layer.partitioning.begin(), layer.partitioning.end(), model.partitioning.begin());

  //    coarse nodes
  model.centroids = layer.coarse_to_fine;
}


void MultiScaleDataMech::build_boundary_nodes()
{
  auto & layer = active_layer();

  // save vertices on interfaces between blocks
  hash_algorithms::ConnectionMap<unordered_set<size_t>> block_face_vertices;
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    const auto & neighbors = face.neighbors();
    const size_t cell1 = neighbors[0];
    const size_t block1 = layer.partitioning[cell1];
    if (neighbors.size() == 1) continue;
    const size_t cell2 = neighbors[1];
    const size_t block2 = layer.partitioning[cell2];;

    if (block1 != block2)
    {
      size_t conn_index;
      if ( block_face_vertices.connection_exists(block1, block2) )
        conn_index = block_face_vertices.connection_index(block1, block2);
      else
        conn_index = block_face_vertices.insert_connection(block1, block2);

      auto & vertices = block_face_vertices.get_data(conn_index);
      for (const size_t vertex : face.vertex_indices())
        vertices.insert(vertex);
    }
  }

  // fill boundary vertices
  // idea: take coarse faces between a block that contains the coarse
  // vertex and its neighbors that do not
  layer.support_boundary.resize(layer.coarse_to_fine.size());
  for (size_t coarse_vertex = 0; coarse_vertex < layer.coarse_to_fine.size(); coarse_vertex++)
  {
    const size_t fine_vertex = layer.coarse_to_fine[coarse_vertex];
    std::cout << "coarse "<<coarse_vertex << " fine " <<fine_vertex << std::endl;
    // auto it = map_coarse_node_blocks.find(fine_vertex);
    // if (it == map_coarse_node_blocks.end())
    // {
    //   std::cout << "SHITTTT" << std::endl;
    //   exit(0);
    // }
    // const auto & neighboring_blocks = it->second;
    const auto & neighboring_blocks = map_coarse_node_blocks[fine_vertex];
    for (const size_t block1 : neighboring_blocks)
    {
      std::cout << "support in block " << block1 << std::endl;
      const auto blocks2 = block_face_vertices.get_neighbors(block1);
      for (const size_t block2 : blocks2)
      {
        std::cout << "\ttrying block " << block2 <<" ";
        if (neighboring_blocks.find(block2) == neighboring_blocks.end())
        {
          std::cout << "OK" << std::endl;
          for (const size_t node : block_face_vertices.get_data(block1, block2))
            layer.support_boundary[coarse_vertex].insert(node);
        }
        else std::cout << "NO" << std::endl;
      }
    }
  }
}

}  // end namespace
