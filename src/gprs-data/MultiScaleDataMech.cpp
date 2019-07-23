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

  algorithms::UnionFindWrapper<size_t> face_disjoint = build_external_face_disjoint();
  size_t ghost_block = active_layer().n_blocks;
  std::unordered_map<size_t, size_t> map_block_group;
  for ( const auto &it_face: face_disjoint.items() )
    if (map_block_group.find(it_face.first) == map_block_group.end())
      map_block_group.insert({ it_face.first, ghost_block++});

  const auto cell_block_neighbors = build_cell_block_neighbors(face_disjoint, map_block_group);
  find_block_corners(face_disjoint, map_block_group, cell_block_neighbors);
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
  unordered_set<size_t> global_corners;
   for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
     if (cell_block_neighbors[cell.index()].size() > 2)  // contains block corner
     {
       // don't want a hash_map here
       map<size_t, unordered_set<size_t>> vertex_blocks;

       // this is by default when we access the element apparently
       // for (const size_t vertex: cell.vertices())
       //   vertex_count.insert( {vertex, 0} );

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
           }
         }
       }

       // there should be only
       for (auto & it : vertex_blocks)
         if (it.second.size() > 2)  // it's a corner!
         {
           if (global_corners.find(it.first) == global_corners.end())
           {
             layer.coarse_node_vertices.push_back(it.first);

           }
           else continue;
         }
     }

}

// void MultiScaleDataMech::fill_output_model(MultiScaleOutputData & model, const int layer_index)
// {

// }

}  // end namespace
