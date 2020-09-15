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
  std::cout << "building METIS partitioning...";
  build_partitioning();
  std::cout << "OK" << std::endl;
  build_cells_in_block();

  const auto map_boundary_face_ghost_block = build_map_face_to_ghost_cell();
  const auto cell_block_neighbors = build_cell_block_neighbors(map_boundary_face_ghost_block);
  find_block_corners2(map_boundary_face_ghost_block, cell_block_neighbors);
  build_boundary_nodes(map_boundary_face_ghost_block);
}


std::vector<std::unordered_set<std::size_t>>
// MultiScaleDataMech::build_cell_block_neighbors(const algorithms::UnionFindWrapper<size_t> & face_disjoint,
//                                                const std::unordered_map<size_t, size_t>   & map_block_group) const
MultiScaleDataMech::
build_cell_block_neighbors(const std::unordered_map<size_t, size_t> & map_boundary_face_ghost_block) const
{
  const auto & layer = active_layer();
  vector<unordered_set<size_t>> cell_block_neighbors(grid.n_cells_total());

  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
  {
    const auto & neighbors = face->neighbors();
    const size_t cell1 = neighbors[0]->index();
    const size_t block1 = layer.partitioning[cell1];
    if (neighbors.size() == 2)
    {
      const size_t cell2 = neighbors[1]->index();
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
      const size_t block2 = map_boundary_face_ghost_block.find(face->index())->second;
      cell_block_neighbors[cell1].insert(block1);
      cell_block_neighbors[cell1].insert(block2);
    }
  }
  return cell_block_neighbors;
}


void MultiScaleDataMech::
find_block_corners2(const std::unordered_map<size_t, size_t> & map_boundary_face_ghost_block,
                    const std::vector<std::unordered_set<std::size_t>> & cell_block_neighbors)
{
  auto & layer = active_layer();

  std::unordered_map<size_t, std::unordered_set<size_t>> map_coarse_node_blocks;
  for (std::size_t block=0; block<layer.n_blocks; ++block)
  {
    unordered_map<size_t, unordered_set<size_t>> vertex_blocks;
    // unordered_map<size_t, unordered_set<size_t>> vertex_faces;

    for (const size_t i : layer.cells_in_block[block])
     if (cell_block_neighbors[i].size() > 1)  // might contain block corner
     {
       const auto & cell = grid.cell(i);
       // const size_t block1 = layer.partitioning[i];
       // assert(block1 == block);
       // a corner is a node all of whose adjacent faces touch different blocks
       for (const auto & face : cell.faces())
       {
         size_t block2;
         const auto & neighbors = face->neighbors();
         if (neighbors.size() == 2)
         {
           if (neighbors[0]->index() == i)
             block2 = layer.partitioning[neighbors[1]->index()];
           else if (neighbors[1]->index() == i)
             block2 = layer.partitioning[neighbors[0]->index()];
           else
             assert(false);
         }
         else  // if (neighbors.size() == 1)  -- ghost case
           block2 = map_boundary_face_ghost_block.find(face->index())->second;

         if (block != block2)
           for (const size_t vertex : face->vertices())
           {
             vertex_blocks[vertex].insert(block);
             vertex_blocks[vertex].insert(block2);
           }

       }  // end cell face loop
     }  // end cell loop within block

    // if more than two vertex blocks
    for (auto & it : vertex_blocks)
      if (it.second.size() > 3)
      { // it's a corner!
        auto glob_it = map_coarse_node_blocks.find(it.first);
        if (glob_it == map_coarse_node_blocks.end())
          map_coarse_node_blocks.insert({ it.first, it.second});
        else
        {
          auto & blocks = glob_it->second;
          for (const size_t b : it.second) blocks.insert(b);
        }
      }
  }

   //  save coarse node vertices
   layer.coarse_to_fine.resize(map_coarse_node_blocks.size());
   coarse_node_blocks.resize(map_coarse_node_blocks.size());
   size_t index = 0;
   for (const auto & it : map_coarse_node_blocks)
   {
     layer.coarse_to_fine[index] = it.first;
     for (const size_t block : it.second)
       // if (!is_ghost_block(block))
         coarse_node_blocks[index].push_back(block);

     index++;
   }
}


void MultiScaleDataMech::
find_block_corners(const std::unordered_map<size_t, size_t> & map_boundary_face_ghost_block,
                   const std::vector<std::unordered_set<std::size_t>> & cell_block_neighbors)
{
  auto & layer = active_layer();
  std::unordered_map<size_t, std::unordered_set<size_t>> map_coarse_node_blocks;
  // given a list of cells that have more than two neighboring blocks
  // find corners of the blocks.
  // To do this we loop over the faces of the cell, figure out which faces touch
  // other blocks, and then find a common vertex among those vertices
  // unordered_set<size_t> global_corners;
   for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
     if (cell_block_neighbors[cell->index()].size() > 2)  // might contain block corner
     {
       // don't want a hash_map here
       map<size_t, unordered_set<size_t>> vertex_blocks;
       map<size_t, unordered_set<size_t>> vertex_faces;

       // a corner is a node all of whose adjacent faces touch different blocks
       for (const auto & face : cell->faces())
       {
         const auto & neighbors = face->neighbors();
         const size_t cell1 = neighbors[0]->index();
         const size_t block1 = layer.partitioning[cell1];
         size_t block2;
         if (neighbors.size() == 2)
           block2 = layer.partitioning[neighbors[1]->index()];
         else // if (neighbors.size() == 1)
           block2 = map_boundary_face_ghost_block.find(face->index())->second;

         if (block1 != block2)
           for (const auto & vertex : face->vertices())
           {
             vertex_blocks[vertex].insert(block2);
             vertex_faces[vertex].insert(face->index());
           }
       }

       for (auto & it : vertex_blocks)
         // technically we should check that ALL vertex faces
         // are touching different blocks.
         // But this criteria is enough for tetras and hexes
         if (it.second.size() > 2)
           if (vertex_faces[it.first].size() > 2)
           { // it's a corner!
             auto glob_it = map_coarse_node_blocks.find(it.first);
             if (glob_it == map_coarse_node_blocks.end())
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

   // save coarse node vertices
   layer.coarse_to_fine.resize(map_coarse_node_blocks.size());
   size_t index = 0;
   for (const auto & it : map_coarse_node_blocks)
   {
     std::cout << "index = " << index << std::endl;
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

  // cupport boundary
  model.support_boundary = layer.support_boundary;

  // support internal nodes are all the block nodes except for boundary nodes
  model.support_internal.resize(layer.coarse_to_fine.size());
  for (size_t coarse_vertex = 0; coarse_vertex < layer.coarse_to_fine.size(); coarse_vertex++)
  {
    const size_t fine_vertex = layer.coarse_to_fine[coarse_vertex];
    const auto & neighboring_blocks = coarse_node_blocks[coarse_vertex];

    // approximate number of nodes to allocate memory
    // size_t n_approx_internal_cells = 0;
    size_t n_approx_internal_cells = 0;
    for (const size_t block : neighboring_blocks)
      if (!is_ghost_block(block))
        n_approx_internal_cells += layer.cells_in_block[block].size();

    // fill internal cells
    model.support_internal[coarse_vertex].reserve(n_approx_internal_cells);
    for (const size_t block : neighboring_blocks)
      if (!is_ghost_block(block))
        for(const size_t cell: layer.cells_in_block[block])
          model.support_internal[coarse_vertex].insert(cell);
  }

  std::cout << std::endl;
  std::cout << "#########################" << std::endl;
  std::cout << "Exported multiscale model" << std::endl;
  std::cout << "Partitioning size: " << layer.n_blocks << std::endl;
  std::cout << "Number of coarse nodes: " << model.n_coarse << std::endl;
  // for (size_t coarse_vertex = 0; coarse_vertex < layer.coarse_to_fine.size(); coarse_vertex++)
    // std::cout << coarse_vertex << " "
    //           << model.centroids[coarse_vertex] << " "
    //           << model.support_internal[coarse_vertex].size() << " "
    //           << model.support_boundary[coarse_vertex].size() << " "
    //           << std::endl;
  std::cout << "#########################" << std::endl;
  std::cout << std::endl;

}


void MultiScaleDataMech::build_boundary_nodes(const std::unordered_map<size_t, size_t> & map_boundary_face_ghost_block)
{
  auto & layer = active_layer();

  // save fine vertices on coarse between blocks
  hash_algorithms::ConnectionMap<unordered_set<size_t>> block_face_vertices;
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
  {
    const auto & neighbors = face->neighbors();
    const size_t cell1 = neighbors[0]->index();
    const size_t block1 = layer.partitioning[cell1];
    size_t block2;
    if (neighbors.size() == 1)
      block2 = map_boundary_face_ghost_block.find(face->index())->second;
    else // if (neighbors.size() == 2)
    {
      const size_t cell2 = neighbors[1]->index();
      block2 = layer.partitioning[cell2];;
    }

    if (block1 != block2)
    {
      size_t conn_index;
      if ( block_face_vertices.contains(block1, block2) )
        conn_index = block_face_vertices.index(block1, block2);
      else
        conn_index = block_face_vertices.insert(block1, block2);

      auto & vertices = block_face_vertices.get_data(conn_index);
      for (const size_t vertex : face->vertices())
        vertices.insert(vertex);
    }
  }

  // fill support region boundary vertices
  // take coarse faces between a block that contains the coarse
  // vertex and its neighbors that do not
  layer.support_boundary.resize(layer.coarse_to_fine.size());

  for (size_t coarse_vertex = 0; coarse_vertex < layer.coarse_to_fine.size(); coarse_vertex++)
  {
    const size_t fine_vertex = layer.coarse_to_fine[coarse_vertex];

    const auto & neighboring_blocks = coarse_node_blocks[coarse_vertex];
    for (const size_t block1 : neighboring_blocks)
      if (!is_ghost_block(block1))
    {
      const auto blocks2 = block_face_vertices.get_neighbors(block1);
      for (const size_t block2 : blocks2)
      {
        if (std::find(neighboring_blocks.begin(), neighboring_blocks.end(), block2) ==
            neighboring_blocks.end())
        {
          for (const size_t node : block_face_vertices.get_data(block1, block2))
            layer.support_boundary[coarse_vertex].insert(node);
        }
      }
    }
  }
}

}  // end namespace
