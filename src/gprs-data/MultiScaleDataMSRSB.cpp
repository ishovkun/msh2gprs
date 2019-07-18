#include "MultiScaleDataMSRSB.hpp"
#include "MetisInterface.hpp"

#include <unordered_set>

namespace multiscale {

using Point = angem::Point<3,double>;
using std::unordered_set;

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
    PureConnectionMap cell_connections;

    for (auto it = grid.begin_faces(); it != grid.end_faces(); ++it)
    {
      const auto & neighbors = it.neighbors();
      if (neighbors.size() == 2)  // not a boundary face
        cell_connections.insert_connection( neighbors[0], neighbors[1] );
    }

    layer.partitioning = multiscale::MetisInterface<hash_algorithms::empty>
        ::build_partitioning(cell_connections, layer.n_blocks, layer.n_cells);
}


void MultiScaleDataMSRSB::build_support_regions()
{
  find_centroids();
  build_block_connections();
  build_block_face_data();
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


void MultiScaleDataMSRSB::build_block_connections()
{
  auto & layer = active_layer();

  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    const auto & neighbors = face.neighbors();
    if (neighbors.size() == 2)  // not a boundary face
    {
      const std::size_t i1 = layer.partitioning[neighbors[0]];
      const std::size_t i2 = layer.partitioning[neighbors[1]];;
      if (i1 != i2)
        if (!layer.block_internal_connections.connection_exists(i1, i2))
          layer.block_internal_connections.insert_connection(i1, i2);
    }
  }
}

void MultiScaleDataMSRSB::build_block_face_data()
{
  algorithms::UnionFindWrapper<size_t> face_disjoint = build_external_face_disjoint();

  size_t ghost_block = active_layer().n_blocks;
  std::unordered_map<size_t, size_t> map_block_group;
  for ( const auto &it_face: face_disjoint.items() )
    if (map_block_group.find(it_face.first) == map_block_group.end())
      map_block_group.insert({ it_face.first, ghost_block++});

  find_block_face_centroids(face_disjoint, map_block_group);
  const auto map_vertex_blocks = build_map_vertex_blocks(face_disjoint, map_block_group);
  const auto map_block_vertices = build_block_corners(map_vertex_blocks);
}


void MultiScaleDataMSRSB::find_block_face_centroids(algorithms::UnionFindWrapper<size_t> & face_disjoint,
                                                    std::unordered_map<size_t, size_t>   & map_block_group)
{
  auto & layer = active_layer();
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    const auto & neighbors = face.neighbors();
    size_t block1, block2;
    block1 = layer.partitioning[neighbors[0]];

    if (neighbors.size() == 1)  // block + ghost block
    {
      const size_t face_group = face_disjoint.group(face.index());
      const size_t block2 = map_block_group[face_group];
    }
    else // if (neighbors.size() == 2)
    {
      block2 = layer.partitioning[neighbors[1]];
      if (block1 == block2) continue;
    }

    size_t conn_ind;
    if (layer.block_faces.connection_exists(block1, block2))
      conn_ind = layer.block_faces.connection_index(block1, block2);
    else
      conn_ind = layer.block_faces.insert_connection(block1, block2);

    auto & data = layer.block_faces.get_data(conn_ind);
    data.center += face.center();
    data.n_cell_faces++;
  }

  for (auto face = layer.block_faces.begin(); face != layer.block_faces.end(); ++face)
    face->center /= face->n_cell_faces;
}


bool MultiScaleDataMSRSB::share_edge(const mesh::const_face_iterator &face1,
                                     const mesh::const_face_iterator &face2)
{
  unordered_set<size_t> shared_verts;
  for (const auto &v1 : face1.vertex_indices())
    for (const auto &v2 : face2.vertex_indices())
      if (v1 == v2)
      shared_verts.insert(v1);

  return (shared_verts.size() > 1);
}


algorithms::UnionFindWrapper<size_t> MultiScaleDataMSRSB::build_external_face_disjoint()
{
  auto & layer = active_layer();

  // build map vertex - external face and fill disjoints
  std::unordered_map<size_t, vector<mesh::const_face_iterator>> map_vertex_face;
  algorithms::UnionFindWrapper<size_t> face_disjoint;

  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
    if (face.neighbors().size() == 1)
    {
      face_disjoint.insert(face.index());

      for (size_t & v : face.vertex_indices() )
      {
        auto it = map_vertex_face.find(v);
        if (it == map_vertex_face.end()) map_vertex_face.insert({v, {face}});
        else it->second.push_back(face);
      }
    }

  // create external face groups
  for (const auto & it : map_vertex_face)
  {
    const auto & faces = it.second;
    for (std::size_t i=0; i<faces.size(); ++i)
    {
      const auto face1 = faces[i];
      const size_t cell1 = face1.neighbors()[0];
      for (std::size_t j=i+1; j<faces.size(); ++j)
      {
        const auto face2 = faces[j];
        const size_t cell2 = face2.neighbors()[0];
        // if (face1 != face2) // they are created different

        if (cell1 != cell2)  // criterion 1 for groups
          if ( abs(dot(face1.normal(), face2.normal())) > 1e-4 ) // criterion 2
            // note: share not just vertex but an edge
            if ( share_edge(face1, face2) )                      // criterion 3
          {
            const size_t block1 = layer.partitioning[cell1];
            const size_t block2 = layer.partitioning[cell2];
            face_disjoint.merge(face1.index(), face2.index());
          }

      }
    }
  }
  return face_disjoint;
}


std::unordered_map<std::size_t, std::vector<std::size_t>>
MultiScaleDataMSRSB::
build_map_vertex_blocks(algorithms::UnionFindWrapper<size_t> & face_disjoint,
                        std::unordered_map<size_t, size_t>   & map_block_group)
{
  auto & layer = active_layer();
  std::unordered_map<std::size_t, std::vector<std::size_t>> vertex_blocks;
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    const auto & neighbors = face.neighbors();
    size_t block1, block2;
    block1 = layer.partitioning[neighbors[0]];

    if (neighbors.size() == 2)
    {
      block2 = layer.partitioning[neighbors[1]];
      if (block1 == block2) continue;
    }
    else if (neighbors.size() == 1)
    {
      const size_t face_group = face_disjoint.group(face.index());
      block2 = map_block_group[face_group];
    }

    for (const auto & vertex : face.vertex_indices())
    {
      auto it_vert = vertex_blocks.find(vertex);

      if (it_vert == vertex_blocks.end())
        vertex_blocks.insert({vertex, {block1, block2}});
      else
      {
        auto & blocks = it_vert->second;
        // std::find works here since a vertex probaby doesn't connect to
        // more than a few blocks
        if (std::find(blocks.begin(), blocks.end(), block1) == blocks.end())
          blocks.push_back(block1);
        if (std::find(blocks.begin(), blocks.end(), block2) == blocks.end())
          blocks.push_back(block2);
      }

    }
  }
  return vertex_blocks;
}

bool MultiScaleDataMSRSB::is_ghost_block(const size_t block) const
{
  if (block < active_layer().n_blocks)
    return false;
  else return true;
}

}  // end namespace
