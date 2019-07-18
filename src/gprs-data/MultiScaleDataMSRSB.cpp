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

  std::cout << "building METIS partitioning...";
  build_partitioning();
  std::cout << "OK" << std::endl;
  build_support_regions();
}


void MultiScaleDataMSRSB::build_partitioning()
{
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
  std::cout << "finding block centroids...";
  find_centroids();
  std::cout << "OK" << std::endl;

  std::cout << "building block connections...";
  build_block_connections();
  std::cout << "OK" << std::endl;

  std::cout << "build block face data..." << std::flush;
  build_block_face_data();
  std::cout << "OK" << std::endl;
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
  std::cout << "building disjoints..." << std::flush;
  algorithms::UnionFindWrapper<size_t> face_disjoint = build_external_face_disjoint();
  std::cout << "OK" << std::endl;

  size_t ghost_block = active_layer().n_blocks;
  std::unordered_map<size_t, size_t> map_block_group;
  for ( const auto &it_face: face_disjoint.items() )
    if (map_block_group.find(it_face.first) == map_block_group.end())
      map_block_group.insert({ it_face.first, ghost_block++});

  std::cout << "building face centroids...";
  find_block_face_centroids(face_disjoint, map_block_group);
  std::cout << "OK" << std::endl;
  std::cout << "building map vert blocks ...";
  const auto map_vertex_blocks = build_map_vertex_blocks(face_disjoint, map_block_group);
  std::cout << "OK" << std::endl;
  std::cout << "finding blocks corners ...";
  const auto map_block_edge_vertices = build_block_edges(map_vertex_blocks);
  std::cout << "OK" << std::endl;
  std::cout << "find block edge centers...";
  find_block_edge_centroids(map_block_edge_vertices);
  std::cout << "OK" << std::endl;
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
        if (it == map_vertex_face.end()) map_vertex_face.insert({v, {std::move( face )}});
        else it->second.push_back(std::move( face ));
      }
    }

  face_disjoint.finalize();

  //  create external face groups
  for (const auto & it : map_vertex_face)
  {
    const auto & faces = it.second;
    for (std::size_t i=0; i<faces.size(); ++i)
    {
      const auto & face1 = faces[i];
      const size_t cell1 = face1.neighbors()[0];
      for (std::size_t j=i+1; j<faces.size(); ++j)
      {
        const auto & face2 = faces[j];
        const size_t cell2 = face2.neighbors()[0];
        // if (face1 != face2) // they are created different

        /* criteria for uniting faces into groups to identify
         * ghost blocks:
         * (1) two boundary face do not belong to the  same cell
         * (2) their normals don't have a crazy angle between 'em
         * (3) two faces share an edge (two connected vertices) */
        if (cell1 != cell2)                                      // criterion 1
          if ( abs(dot(face1.normal(), face2.normal())) > 1e-4 ) // criterion 2
            if ( share_edge(face1, face2) )                      // criterion 3
            face_disjoint.merge(face1.index(), face2.index());
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
  /* collect vertices from faces that are on block-block interfaces */
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


std::unordered_map<std::tuple<std::size_t,std::size_t,std::size_t>, std::vector<std::size_t>>
MultiScaleDataMSRSB::
build_block_edges(const std::unordered_map<std::size_t, std::vector<std::size_t>> & map_block_vertices)
{
  std::unordered_map<std::tuple<std::size_t,std::size_t,std::size_t>,
                     std::vector<std::size_t>> block_triplet_vertices;

  for (auto & it_vertex : map_block_vertices)
    if (it_vertex.second.size() > 2)  // on edge
    {
      // copy cause modifying
      auto blocks = it_vertex.second;
      // build block triplets with ascending ordering
      std::sort(blocks.begin(), blocks.end());

      for (std::size_t i = 0; i < blocks.size(); ++i)
        for (std::size_t j = i+1; j < blocks.size(); ++j)
          for (std::size_t k = j+1; k < blocks.size(); ++k)
          {
            const auto block_triplet = std::make_tuple(blocks[i], blocks[j], blocks[k]);
            auto it_blocks = block_triplet_vertices.find(block_triplet);
            if (it_blocks != block_triplet_vertices.end())
              it_blocks->second.push_back(it_vertex.first);
            else
              block_triplet_vertices.insert({std::move( block_triplet ), {it_vertex.first}});
          }
    }

  return block_triplet_vertices;
}


void MultiScaleDataMSRSB::find_block_edge_centroids(
      const std::unordered_map<std::tuple<std::size_t,std::size_t,std::size_t>, std::vector<std::size_t>>
      &map_block_edge_vertices)
{
  for (const auto & block_edge : map_block_edge_vertices)
  {
    // find block edge center
    angem::Point<3,double> block_edge_center = {0, 0, 0};
    for (const auto & vertex : block_edge.second)
      block_edge_center += grid.vertex_coordinates(vertex);
    block_edge_center /= block_edge.second.size();

    // move this point to lie on the edge
    // make this point lie on edge
    // block_edge_center = grid.vertex_coordinates(
    //     angem::find_closest_vertex(block_edge_center,
    //           /* all_vertices = */ grid.get_vertices(),
    //                 /* subset = */ block_edge.second));

    // std::pair<std::size_t, std::size_t> face_blocks;
    // std::size_t edging_block;
  }
}

}  // end namespace
