#include <Mesh.hpp>
#include <SurfaceMesh.hpp>
#include "Cell.hpp"
#include "angem/PolyhedronFactory.hpp"
#include "angem/Collisions.hpp"    // angem::split
#include <unordered_set>
#include <algorithm>  // std::max

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
  return insert_cell_(ivertices, poly_faces, vtk_id, marker);
}

std::size_t Mesh::
insert_cell_(const std::vector<std::vector<std::size_t>> & cell_faces,
            const int                        marker)
{
  // cells don't have many vertices: don't need a hashmap here
  std::set<size_t> unique_vertices;
  for (const auto & face : cell_faces)
    for (const size_t v : face)
      unique_vertices.insert(v);
  std::vector<size_t> ivertices(unique_vertices.begin(), unique_vertices.end());
  return insert_cell_( ivertices, cell_faces, -1, marker );
}

std::size_t Mesh::
insert_cell_(const std::vector<std::size_t> & ivertices,
             const std::vector<std::vector<std::size_t>> & cell_faces,
             const int                        vtk_id,
             const int                        marker)
{
  const std::size_t new_cell_index = n_cells();

  std::vector<std::size_t> face_indices;
  face_indices.reserve(cell_faces.size());

  for (const vector<size_t> & face : cell_faces)
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

  for (const size_t vertex : ivertices)
    if (std::find( m_vertex_faces[vertex].begin(), m_vertex_faces[vertex].end(),
                   face_index) == m_vertex_faces[vertex].end())
      m_vertex_faces[vertex].push_back(face_index);

  return face_index;
}


void Mesh::split_vertex(const std::size_t               vertex_index,
                        const std::vector<std::size_t> &splitted_face_indices)
{
  // this code is for when splitted face indices contains all different faces
  // I'm not really doing anything with the vertices,
  // cause adgprs doesn't need that
  std::vector<std::vector<std::size_t>> groups =
      group_cells_based_on_split_faces(m_vertex_cells[vertex_index],
                                       splitted_face_indices);
  // create new vertices
  std::vector<std::size_t> new_vertex_indices(groups.size());
  const angem::Point<3,double> vertex_coord = m_vertices[vertex_index];
  for (std::size_t group = 0; group < groups.size(); group++)
  {
    if (group == 0)  // group 0 retains old vertex
      new_vertex_indices[group] = vertex_index;
    else  // add new vertices
    {
      const std::size_t new_vertex_index = m_vertices.size();
      new_vertex_indices[group] = new_vertex_index;
      m_vertices.push_back(vertex_coord);
    }
  }

  // modify cell vertices: replace vertex indices with the new vertices
  // start from 1 since group 0 retains the old index
  for (std::size_t group = 1; group < groups.size(); group++)
  {
    const std::vector<size_t> & cell_group = groups[group];
    for (const std::size_t cell_index : cell_group)
    {
      std::vector<size_t> & cell_vertices = m_cells[cell_index].vertices();
      for (size_t & cell_vertex_index : cell_vertices)
        if (cell_vertex_index == vertex_index)
        {
          cell_vertex_index = new_vertex_indices[group];
        }

    }
  }


}


void Mesh::mark_for_split(const std::size_t face_index)
{
  assert( face_index < n_faces() );
  m_faces_marked_for_split.push_back(face_index);
}


SurfaceMesh<double> Mesh::split_faces()
{
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
    split_vertex(vertex, faces);
  }

  return mesh_faces;
}


std::vector<std::vector<std::size_t>>
Mesh::
group_cells_based_on_split_faces(const std::vector<size_t> & affected_cells,
                                 const std::vector<size_t> & split_faces) const
{
  // group affected elements
  // two elements are in the same group if they are neighbors and
  // the neighboring face is not in vertex_faces array
  const size_t n_groups = std::max(split_faces.size(), size_t(2));
  std::unordered_map<std::size_t, size_t> map_cell_group;
  int igroup = 0;
  int new_group = 0;
  std::unordered_set<std::size_t> processed_cells;
  for (const std::size_t icell : affected_cells)
  {
    auto group_it = map_cell_group.find(icell);
    if (group_it != map_cell_group.end())
      igroup = group_it->second;
    else
    {
      igroup = new_group;
      map_cell_group.insert({icell, igroup});
      new_group++;
    }

    processed_cells.insert(icell);
    // std::cout << "\nicell = " << icell << std::endl;

    // find neighboring cell from affected cells group
    for (const Cell* jcell : m_cells[icell].neighbors())
    {
      if (std::find(affected_cells.begin(), affected_cells.end(),
                    jcell->index()) != affected_cells.end())
      {
        // take index explicitly since minmax takes a reference
        const size_t jind = jcell->index();
        // what face neighbors should be
        const auto pair_cells = std::minmax(icell, jind);

        // find out if i and j neighbor by a marked face
        bool neighbor_by_marked_face = false;
        for (const size_t iface : split_faces)
        {
          const Face f = face(iface);
          const auto f_neighbors = f.neighbors();
          assert( f_neighbors.size() == 2 );
          auto pair_cells2 = std::minmax(f_neighbors[0]->index(),
                                         f_neighbors[1]->index());

          if (pair_cells == pair_cells2)
          // if (pair_cells.first == pair_cells2.first &&
          //     pair_cells.second == pair_cells2.second)
          {
            neighbor_by_marked_face = true;
            break;
          }
        }

        if (!neighbor_by_marked_face)
        {
          auto group_it = map_cell_group.find(jcell->index());
          if (group_it == map_cell_group.end())
            map_cell_group.insert({jcell->index(), igroup});
          else
          {
            if (group_it->second < igroup)
            {
              map_cell_group[icell] = group_it->second;
              igroup = group_it->second;
              new_group--;
            }
            else
              map_cell_group[jcell->index()] = igroup;
          }
        }
      }
    }
  }

  std::vector<std::vector<std::size_t>> groups(n_groups);
  for (auto it : map_cell_group)
    groups[it.second].push_back(it.first);

  return groups;
}


void Mesh::split_cell(Cell & cell,
                      const angem::Plane<double> & plane)
{
  const std::unique_ptr<angem::Polyhedron<double>> polyhedron = cell.polyhedron();
  // Bookkeeping:
  //  fill polygroup's internal set with the existing vertex coordinates
  // in order to have a map of those to the global vertex indices,
  // which will come in handy when inserting new splitted cells into grid.
  // We can do it because splitting will insert the same vertices plus
  // those that appeared due to plase-face intersection.
  angem::PolyGroup<double> split;
  std::vector<size_t> global_vertex_indices;
  for (const Point & p : polyhedron->get_points())
    global_vertex_indices.push_back(split.vertices.insert(p));

  // the actual geometry happens here
  angem::split(*polyhedron, plane, split,
               constants::marker_below_splitting_plane,
               constants::marker_above_splitting_plane,
               constants::marker_splitting_plane);

  // check we actually split something
  assert( split.vertices.size() > global_vertex_indices.size() );

  // insert new vertices (those that occured due to splitting)
  for (size_t i = global_vertex_indices.size(); i < split.vertices.size(); ++i)
  {
    const size_t new_vertex_index = n_vertices();
    m_vertices.push_back(split.vertices[i]);
    global_vertex_indices.push_back(new_vertex_index);
  }

  // make two groups of faces (polyhedra) that will form the new cells
  vector<vector<size_t>> cell_above_faces, cell_below_faces;
  for (size_t i = 0; i < split.polygons.size(); i++)
  {
    if ( split.markers[i] == constants::marker_below_splitting_plane ||
         split.markers[i] == constants::marker_splitting_plane )
      cell_below_faces.push_back(split.polygons[i]);
    if ( split.markers[i] == constants::marker_above_splitting_plane ||
         split.markers[i] == constants::marker_splitting_plane )
      cell_above_faces.push_back(split.polygons[i]);
  }

  // insert new cells
  const size_t child_cell_index1 = insert_cell_(cell_above_faces, cell.marker());
  const size_t child_cell_index2 = insert_cell_(cell_below_faces, cell.marker());

  // handle parent/child dependencies
  cell.m_children = {child_cell_index1, child_cell_index2};
  m_cells[child_cell_index1].m_parent = cell.index();
  m_cells[child_cell_index2].m_parent = cell.index();
}

active_cell_iterator Mesh::begin_active_cells()
{
  for (auto cell = begin_cells(); cell != end_cells(); ++ cell)
    if (!cell->is_active())
      cell++;
    else
      return active_cell_iterator(&*cell);
  return active_cell_iterator(nullptr);
}

}  // end namespace mesh
