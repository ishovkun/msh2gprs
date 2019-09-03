#include <Mesh.hpp>
#include <SurfaceMesh.hpp>
#include "Cell.hpp"
#include "angem/PolyhedronFactory.hpp"
#include <unordered_set>

namespace mesh
{

// using std::unordered_map;
using std::vector;

Mesh::Mesh()
{}


std::size_t Mesh::insert_cell(const std::vector<std::size_t> & ivertices,
                              const int                        vtk_id,
                              const int                        marker)
{
  // get faces: vectors of vertexindices based on vtk numbering
  const std::vector<std::vector<std::size_t>> poly_faces =
      angem::PolyhedronFactory::get_global_faces<double>(ivertices, vtk_id);

  const std::size_t new_cell_index = n_cells();

  std::vector<std::size_t> face_indices;
  face_indices.reserve(poly_faces.size());

  for (const vector<size_t> & face : poly_faces)
  {
    size_t face_vtk_id;
    switch (face.size())
    {
      case 3:
        face_vtk_id = 5;  //  triangle
        break;
      case 4:
        face_vtk_id = 9;  //  vtk_quad
        break;
      default:
        face_vtk_id = 7;  //  vtk_polygon
        break;
    }

    const std::size_t face_index = insert_face(face, face_vtk_id, DEFAULT_FACE_MARKER);
    face_indices.push_back(face_index);
  }

  m_cells.emplace_back(new_cell_index, ivertices, std::move(face_indices),
                       m_vertices, m_cells, m_faces, vtk_id, marker);

  for (const size_t vertex: ivertices)
    m_vertex_cells[vertex].push_back(new_cell_index);

  return new_cell_index;
}


size_t Mesh::insert_face(const std::vector<std::size_t> & ivertices,
                         const int                        vtk_id,
                         const int                        marker)
{
  if (m_vertex_cells.size() < n_vertices())
    m_vertex_cells.resize(n_vertices());
  if (m_vertex_faces.size() < n_vertices())
    m_vertex_faces.resize(n_vertices());

  size_t face_index = m_faces.size();
  std::vector<size_t> sorted_vertices(ivertices);
  std::sort(sorted_vertices.begin(), sorted_vertices.end());

  bool face_exists = false;
  for (const size_t vertex: ivertices)
  {
    for (const size_t iface: m_vertex_faces[vertex])
    {
      assert( n_faces() > iface );
      Face & f = m_faces[iface];
      std::vector<size_t> face_vertices(f.vertices());
      std::sort(face_vertices.begin(), face_vertices.end());
      if ( std::equal( face_vertices.begin(), face_vertices.end(),
                       sorted_vertices.begin(), sorted_vertices.end()) )
      {
        face_exists = true;
        face_index = iface;
        break;
      }
    }
  }

  if (!face_exists)
  {
    m_faces.emplace_back(face_index, face_index,
                         ivertices, vtk_id, marker,
                         m_cells, m_vertices, m_vertex_cells);
  }
  else
  {
    if (marker != DEFAULT_FACE_MARKER)
      m_faces[face_index].m_marker = marker;
    assert( m_faces[face_index].m_vtk_id == vtk_id );
  }

  if (face_index == 0 ||
      (  sorted_vertices[0] == 10 and
         sorted_vertices[1] == 22 and
         sorted_vertices[2] == 270))
  {
    std::cout << "\n "<< "##########" << std::endl;
    std::cout << "face_index = " << face_index << std::endl;
      std::cout << "vertices:" << std::endl;
      for (auto v: ivertices)
        std::cout << v << " ";
      std::cout << std::endl;
  }

  for (const size_t vertex : ivertices)
    if (std::find( m_vertex_faces[vertex].begin(), m_vertex_faces[vertex].end(),
                   face_index) == m_vertex_faces[vertex].end())
      m_vertex_faces[vertex].push_back(face_index);

  return face_index;
}


void Mesh::split_vertex(const std::size_t vertex_index,
                        const std::size_t master_face_index)
{
  std::vector<std::size_t> affected_face_indices;
  for (const size_t iface : m_vertex_faces[vertex_index])
  {
    const Face & f = m_faces[iface];
    if (f.master_index() == master_face_index)
      if (f.has_vertex(vertex_index))
        affected_face_indices.push_back(f.index());
  }
  
  // for (auto id : affected_face_indices)
  //   std::cout << id << std::endl;
  // exit(0);
  std::vector<std::size_t> affected_cell_indices;
  affected_cell_indices.reserve(2);
  for (const size_t iface : affected_face_indices)
  {
    const Face & f = m_faces[iface];
    for (const Cell * cell : f.neighbors())
      if (std::find( affected_cell_indices.begin(), affected_cell_indices.end(),
                     cell->index() == affected_cell_indices.end()) )
        affected_cell_indices.push_back(cell->index());
      }
}


void Mesh::mark_for_split(const std::size_t face_index)
{
  assert( face_index < n_faces() );
  m_faces_marked_for_split.push_back(face_index);
}


SurfaceMesh<double> Mesh::split_faces()
{
  std::cout << "n_nodes() = " << n_vertices() << std::endl;
  std::cout << "n_faces = " << n_faces() << std::endl;

  /* Algorithm:
  * create SurfaceMesh from marked faces in order to identify
  * vertices to split those whose edge have >1 neighbors)
  * cross-match vertices in 3d Mesh and Surface mesh. */

  // create surfacemesh and map vertices
  SurfaceMesh<double> mesh_faces(1e-6);
  // map 2d-element -> 3d face hash
  std::unordered_map<std::size_t, std::size_t> map_face_surface_element;
  // map surfacemesh vertex -> 3d mesh vertex
  std::unordered_map<std::size_t, std::size_t> map_frac_vertex_vertex;
  for (const std::size_t face_index: m_faces_marked_for_split)
  {
    const Face & f = face(face_index);
    const std::size_t ielement = mesh_faces.insert(f.polygon());
    map_face_surface_element.insert({ielement, face_index});
    const auto frac_poly = mesh_faces.create_poly_iterator(ielement);
    size_t iv = 0;
    for (const Point & v : f.vertex_coordinates())
    {
      size_t ifv = 0;
      for (const Point & frac_vertex : frac_poly.vertex_coordinates())
      {
        if (v == frac_vertex)
        {
          const size_t iv_global = f.vertices()[iv];
          const size_t ifv_global = frac_poly.vertices()[ifv];
          if (map_frac_vertex_vertex.find(ifv_global) == map_frac_vertex_vertex.end() )
            map_frac_vertex_vertex.insert({ ifv_global, iv_global });
          else
            assert( map_frac_vertex_vertex.find(ifv_global)->second == iv_global );
        }
        ifv++;
      }
      iv++;
    }
  }

  // find vertices to split as those those edges have more than one neighbor
  std::unordered_map<size_t, std::vector<size_t>> vertices_to_split;
  for (auto edge = mesh_faces.begin_edges(); edge !=mesh_faces.end_edges(); ++edge)
  {
    const vector<size_t> & edge_neighbors = edge.neighbors();
    std::vector<size_t> grid_face_indices;
    grid_face_indices.reserve(edge_neighbors.size());
    for (size_t ielement: edge_neighbors)
      grid_face_indices.push_back(map_face_surface_element[ielement]);

    if (edge_neighbors.size() > 1)  // internal edge vertex
    {
      const auto edge_vertices = edge.vertex_indices();
      const size_t v1 = map_frac_vertex_vertex[edge_vertices.first];
      const size_t v2 = map_frac_vertex_vertex[edge_vertices.second];
      auto it1 = vertices_to_split.find(v1);
      if (it1 == vertices_to_split.end())
        vertices_to_split.insert({v1, edge_neighbors});
      else
        for (const size_t face : grid_face_indices)
          if ( std::find(it1->second.begin(), it1->second.end(), face) ==
               it1->second.end())
            it1->second.push_back(face);

      auto it2 = vertices_to_split.find(v2);
      if (it2 == vertices_to_split.end())
        vertices_to_split.insert({v2, edge_neighbors});
      else
        for (const size_t face : grid_face_indices)
          if ( std::find(it2->second.begin(), it2->second.end(), face) ==
               it2->second.end())
            it2->second.push_back(face);
    }
  }

  // split 'em
  for (const auto & it : vertices_to_split)
  {
    const size_t vertex = it.first;
    const auto & faces = it.second;
    std::cout << "splitting vertex = " << vertex << std::endl;
    for (const size_t face : faces)
      split_vertex(vertex, face);
  }
  exit(0);
  return mesh_faces;
}


// std::vector<std::vector<std::size_t>>
// Mesh::group_cells_based_on_split_faces(const std::unordered_set<std::size_t> & affected_cells,
//                                        const std::vector<face_iterator>      & vertex_faces) const
// {
//   // group affected elements
//   // two elements are in the same group if they are neighbors and
//   // the neighboring face is not in vertex_faces array
//   const int n_groups = vertex_faces.size();
//   std::unordered_map<std::size_t, int> map_cell_group;
//   int igroup = 0;
//   int new_group = 0;
//   // just is purely for faster checking: cells that are already processed
//   std::unordered_set<std::size_t> processed_cells;
//   for (const std::size_t icell : affected_cells)
//   {
//     auto group_it = map_cell_group.find(icell);
//     if (group_it != map_cell_group.end())
//     {
//       igroup = group_it->second;
//     }
//     else
//     {
//       igroup = new_group;
//       map_cell_group.insert({icell, igroup});
//       new_group++;
//     }

//     processed_cells.insert(icell);

//     for (const std::size_t jcell : get_neighbors(icell))
//       if (affected_cells.find(jcell) != affected_cells.end())
//       {
//         auto pair_cells = std::minmax(icell, jcell);
//         std::vector<std::size_t> ordered_neighbors =
//             {pair_cells.first, pair_cells.second};

//         // find out if i and j neighbor by a marked face
//         bool neighbor_by_marked_face = false;
//         for (const auto & face : vertex_faces)
//         {
//           const auto it = map_faces.find(face.hash());

//           if (it->second.neighbors == ordered_neighbors)
//           {
//             neighbor_by_marked_face = true;
//             break;
//           }
//         }

//         if (!neighbor_by_marked_face)
//         {
//           auto group_it = map_cell_group.find(jcell);
//           if (group_it == map_cell_group.end())
//             map_cell_group.insert({jcell, igroup});
//           else
//           {
//             if (group_it->second < igroup)
//             {
//               map_cell_group[icell] = group_it->second;
//               igroup = group_it->second;
//               new_group--;
//             }
//             else
//               map_cell_group[jcell] = igroup;
//           }
//         }
//         // else
//         //   std::cout << " other" << std::endl;
//       }
//     // std::cout << "igroup = " << igroup << std::endl << std::endl;
//   }
//   // std::cout << "groups" << std::endl;
//   // for (auto & it : map_cell_group)
//   //   std::cout << it.first << "\t" << it.second << std::endl;

//   // abort();

//   std::vector<std::vector<std::size_t>> groups(n_groups);
//   for (auto it : map_cell_group)
//     groups[it.second].push_back(it.first);

//   return groups;
// }


// std::vector<face_iterator> Mesh::get_ordered_faces()
// {
//   std::vector<face_iterator> ordered_faces(this->n_faces(), this->begin_faces());
//   for (auto face=begin_faces(); face!=end_faces(); ++face)
//     ordered_faces[face.index()] = face;
//   return ordered_faces;
// }


//  const angem::Point<3,double> & Mesh::vertex_coordinates(const std::size_t i) const
// {
//   assert(i < n_vertices());
//   return vertices.points[i];
// }


// angem::Point<3,double> & Mesh::vertex_coordinates(const std::size_t i)
// {
//   assert(i < n_vertices());
//   return vertices.points[i];
// }


// void Mesh::add_vertex_to_split(const std::size_t                                        vertex,
//                                const std::vector<std::size_t>                         & edge_neighbors,
//                                const std::unordered_map<std::size_t, size_t>       & map_2d_3d,
//                                std::unordered_set<size_t> & vertices_to_split)
// {
//   auto & faces = vertices_to_split[vertex];
//   for (const auto & neighbor : edge_neighbors)
//   {
//     auto it = create_face_iterator(map_faces.find(map_2d_3d.find(neighbor)->second));
//     if (find(faces.begin(), faces.end(), it) == faces.end())
//       faces.push_back(std::move(it));
//   }

// }


// std::vector<std::size_t> & Mesh::get_vertices(const std::size_t cell)
// {
//   return cells[cell];
// }


// const std::vector<std::size_t> & Mesh::get_vertices(const std::size_t cell) const
// {
//   return cells[cell];
// }


}
