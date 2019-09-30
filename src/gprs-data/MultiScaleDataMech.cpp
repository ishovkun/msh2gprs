#include "MultiScaleDataMech.hpp"
#include <unordered_set>
#include <unordered_map>
#include <map>

namespace multiscale
{

using std::unordered_map;
using std::unordered_set;
using std::map;

MultiScaleDataMech::MultiScaleDataMech(mesh::Mesh  & grid,
                                       const std::array<size_t,3> &  n_blocks,
                                       const PartitioningMethod method,
                                       const size_t elimination_level)
    :
    MultiScaleDataMSRSB(grid, n_blocks, method),
    elimination_level(elimination_level)
{}


void MultiScaleDataMech::build_data()
{
  std::cout << "building partitioning...";
  if (partitioning_method == PartitioningMethod::metis)
    build_metis_partitioning();
  else if (partitioning_method == geometric)
    build_geometric_partitioning();
  std::cout << "OK" << std::endl;

  build_cells_in_block();
  const auto map_boundary_face_ghost_block = build_map_face_to_ghost_cell();
  const auto cell_block_neighbors = build_cell_block_neighbors(map_boundary_face_ghost_block);
  find_block_corners2(map_boundary_face_ghost_block, cell_block_neighbors);
  merge_close_blocks();
  build_boundary_nodes(map_boundary_face_ghost_block);
}




std::vector<std::unordered_set<std::size_t>>
MultiScaleDataMech::
build_cell_block_neighbors(const std::unordered_map<size_t, size_t> & map_boundary_face_ghost_block) const
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
      const size_t block2 = map_boundary_face_ghost_block.find(face.index())->second;
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
       const auto cell = grid.create_const_cell_iterator(i);
       // const size_t block1 = layer.partitioning[i];
       // assert(block1 == block);
       // a corner is a node all of whose adjacent faces touch different blocks
       for (const auto & face : cell.faces())
       {
         size_t block2;
         const auto & neighbors = face.neighbors();
         if (neighbors.size() == 2)
         {
           if (neighbors[0] == i)
           {

             block2 = layer.partitioning[neighbors[1]];
           }
           else if (neighbors[1] == i)
           {
             block2 = layer.partitioning[neighbors[0]];
           }
           else
             assert(false);
         }
         else  // if (neighbors.size() == 1)  -- ghost case
           block2 = map_boundary_face_ghost_block.find(face.index())->second;

         if (block != block2)
           for (const size_t vertex : face.vertex_indices())
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
   for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
     if (cell_block_neighbors[cell.index()].size() > 2)  // might contain block corner
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
         size_t block2;
         if (neighbors.size() == 2)
           block2 = layer.partitioning[neighbors[1]];
         else // if (neighbors.size() == 1)
           block2 = map_boundary_face_ghost_block.find(face.index())->second;

         if (block1 != block2)
           for (const auto & vertex : face.vertex_indices())
           {
             vertex_blocks[vertex].insert(block2);
             vertex_faces[vertex].insert(face.index());
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

  // coarse nodes
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
  std::cout << "#########################" << std::endl;
  std::cout << std::endl;

}


void MultiScaleDataMech::build_boundary_nodes(const std::unordered_map<size_t, size_t> & map_boundary_face_ghost_block)
{
  auto & layer = active_layer();

  // save fine vertices on coarse between blocks
  hash_algorithms::ConnectionMap<unordered_set<size_t>> block_face_vertices;
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    const auto & neighbors = face.neighbors();
    const size_t cell1 = neighbors[0];
    const size_t block1 = layer.partitioning[cell1];
    size_t block2;
    if (neighbors.size() == 1)
      block2 = map_boundary_face_ghost_block.find(face.index())->second;
    else // if (neighbors.size() == 2)
    {
      const size_t cell2 = neighbors[1];
      block2 = layer.partitioning[cell2];;
    }

    if (block1 != block2)
    {
      size_t conn_index;
      if ( block_face_vertices.find(block1, block2) != block_face_vertices.end() )
        conn_index = block_face_vertices.index(block1, block2);
      else
        conn_index = block_face_vertices.insert(block1, block2);

      auto & vertices = block_face_vertices.get_data(conn_index);
      for (const size_t vertex : face.vertex_indices())
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


void MultiScaleDataMech::merge_close_blocks()
{
  auto & layer = active_layer();
  // find tolerance distance
  double cell_size;
  {
    const auto cell0 = grid.begin_cells();
    const auto & verts = cell0.vertices();
    cell_size = verts[1] - verts[0];
  }

  // build node_node neighboring map
  std::vector<std::vector<size_t>> block_nodes(layer.n_blocks + n_ghost);
  for (std::size_t node=0; node<layer.coarse_to_fine.size(); ++node)
    for (size_t block : coarse_node_blocks[node])
      block_nodes[block].push_back(node);

  const auto coarse_node_cells = find_node_neighboring_cells(elimination_level);

  std::unordered_set<size_t> nodes_to_delete;
  for (std::size_t block=0; block<layer.n_blocks; ++block)
  {
    for (std::size_t i=0; i<block_nodes[block].size(); ++i)
    {
      const size_t inode = block_nodes[block][i];
      if(nodes_to_delete.find(inode) == nodes_to_delete.end())
        for (std::size_t j=i+1; j < block_nodes[block].size(); ++j)
        {
          const size_t jnode = block_nodes[block][j];
          if(nodes_to_delete.find(jnode) == nodes_to_delete.end())
          {
            // const angem::Point<3,double> & v1 = grid.vertex( layer.coarse_to_fine[inode] );
            // const angem::Point<3,double> & v2 = grid.vertex( layer.coarse_to_fine[jnode] );
            // if (v1.distance(v2) < 1*cell_size)
            for (const auto & it_i : coarse_node_cells[inode])
            {
              auto it_j = coarse_node_cells[jnode].find(it_i.first);
              if (it_j != coarse_node_cells[jnode].end())
                if (it_j->second == 0)
              {
                assert( inode != jnode );
                nodes_to_delete.insert(jnode);
                const auto & old_blocks = coarse_node_blocks[jnode];
                auto & new_blocks = coarse_node_blocks[inode];
                for (size_t b : old_blocks)
                  if (std::find(new_blocks.begin(), new_blocks.end(), b) == new_blocks.end())
                    new_blocks.push_back(b);
                break;
              }
            }
          }
        }
    }
  }

  // remove from coarse_to_fine and from coarse_node_blocks
  std::vector<std::vector<size_t>> cp_cnb(coarse_node_blocks);
  std::vector<std::size_t> cp_ctf(layer.coarse_to_fine);
  coarse_node_blocks.resize(cp_cnb.size() - nodes_to_delete.size());
  layer.coarse_to_fine.resize(cp_ctf.size() - nodes_to_delete.size());

  size_t ind = 0;
  for (size_t coarse = 0 ; coarse < cp_cnb.size(); coarse++)
    if (nodes_to_delete.find(coarse) == nodes_to_delete.end())
    {
      coarse_node_blocks[ind] = cp_cnb[coarse];
      layer.coarse_to_fine[ind] = cp_ctf[coarse];
      ind++;
    }
  std::cout << "Deleted " << cp_cnb.size() - coarse_node_blocks.size()
            << " nodes that were too close" << std::endl;
}


std::vector<std::unordered_map<size_t,size_t>> MultiScaleDataMech::find_node_neighboring_cells(const size_t level) const
{
  const auto & layer = active_layer();
  std::unordered_map<size_t,size_t> vertex_to_node;
  for (size_t coarse = 0; coarse < layer.coarse_to_fine.size(); coarse++)
    vertex_to_node[layer.coarse_to_fine[coarse]] = coarse;

  // find level 0 neighbors (cells adjacent to node)
  std::vector<std::unordered_map<size_t,size_t>> coarse_node_cells(layer.coarse_to_fine.size());
  for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
    for (const size_t vertex : cell.vertices())
    {
      auto it = vertex_to_node.find(vertex);
      if (it != vertex_to_node.end())
        coarse_node_cells[it->second].insert( {cell.index(), 0} );
    }

  // higher level neighbors: neighbors of neighbors
  for (size_t l = 0; l < level; l++)
    for (std::size_t node=0; node<layer.coarse_to_fine.size(); ++node)
      for (const auto & it : coarse_node_cells[node])
        if (it.second == l)  // take only previous level
        {
          const size_t cell_index = it.first;
          const auto cell = grid.create_const_cell_iterator(cell_index);
          for (const size_t neighbor: cell.neighbor_indices())
            if (coarse_node_cells[node].find(neighbor) == coarse_node_cells[node].end())
              coarse_node_cells[node].insert({ neighbor, l+1 });
        }

  return coarse_node_cells;
}

}  // end namespace
