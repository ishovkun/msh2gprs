#include <Mesh.hpp>
#include <SurfaceMesh.hpp>
#include "Cell.hpp"
#include "angem/PolyhedronFactory.hpp"
#include "angem/Collisions.hpp"    // angem::split
#include <unordered_set>
#include <algorithm>  // std::max
#include  <numeric>   // iota

namespace mesh
{

using std::vector;

Mesh::Mesh()
    : m_n_split_cells(0)
{}

std::size_t Mesh::insert_cell(const std::vector<std::size_t> & ivertices,
                              const int                        vtk_id,
                              const int                        marker)
{
  // get faces: vectors of vertexindices based on vtk numbering
  const std::vector<std::vector<std::size_t>> poly_faces =
      angem::PolyhedronFactory::get_global_faces<double>(ivertices, vtk_id);
  std::vector<FaceTmpData> face_vector(poly_faces.size());
  for (size_t i = 0; i < poly_faces.size(); ++i)
  {
    face_vector[i].vertices = poly_faces[i];
    face_vector[i].vtk_id = face_vtk_id_(poly_faces[i].size());
  }

  std::vector<size_t> take_faces(face_vector.size());
  std::iota(take_faces.begin(), take_faces.end(), 0);
  return insert_cell_(ivertices, take_faces, face_vector,vtk_id, marker);
  // return insert_cell_(ivertices, poly_faces, vtk_id, marker);
}

std::size_t Mesh::
insert_cell_(const std::vector<size_t> take_faces,
             const std::vector<FaceTmpData> &big_face_vector,
             const int                        marker)
{
  // cells don't have many vertices: don't need a hashmap here
  std::set<size_t> unique_vertices;
  // for (const auto & face : cell_faces)
  for (const size_t iface : take_faces)
    // for (const size_t v : big_face_vector)
    for (const size_t v : big_face_vector[iface].vertices)
      unique_vertices.insert(v);

  std::vector<size_t> ivertices(unique_vertices.begin(), unique_vertices.end());
  return insert_cell_(ivertices, take_faces, big_face_vector, constants::vtk_index_general_polyhedron, marker);
}

std::size_t Mesh::
insert_cell_(const std::vector<std::size_t> & ivertices,
             const std::vector<size_t> take_faces,
             const std::vector<FaceTmpData> &big_face_vector,
             const int                        vtk_id,
             const int                        marker)
{
  // parent so that function can be reused for adding cells after splitting
  const std::size_t new_cell_index = n_cells();
  std::vector<std::size_t> face_indices;
  face_indices.reserve(take_faces.size());

  // for (std::size_t i=0; i<take_faces.size(); ++i)
  for (const size_t iface : take_faces)
  {
    const auto & face = big_face_vector[iface];
    size_t face_vtk_id = face.vtk_id;
    const std::size_t face_index = insert_face_(face);
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
                          const int                        marker,
                          const std::size_t                face_parent)
{
  FaceTmpData f;
  f.vertices = ivertices;
  f.vtk_id = vtk_id;
  f.marker = marker;
  f.parent = face_parent;
  return insert_face_(f);
}

size_t Mesh::insert_face_(const FaceTmpData & f)
{
  if (m_vertex_cells.size() < n_vertices())
    m_vertex_cells.resize(n_vertices());
  if (m_vertex_faces.size() < n_vertices())
    m_vertex_faces.resize(n_vertices());

  size_t face_index = m_faces.size();
  std::vector<size_t> sorted_vertices = sort_copy_(f.vertices);

  bool face_exists = false;
  for (const size_t vertex: f.vertices)
  {
    for (const size_t iface: m_vertex_faces[vertex])
    {
      assert( n_faces() > iface );
      const Face & f = m_faces[iface];
      std::vector<size_t> face_vertices = sort_copy_(f.vertices());
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
    m_faces.emplace_back(face_index, f.vertices, f.vtk_id, f.marker,
                         m_cells, m_faces, m_vertices, m_vertex_cells, f.parent);
    if (f.parent != constants::invalid_index)
        m_faces[f.parent].m_children.push_back(face_index);
  }
  else
  {
    if (f.marker != constants::default_face_marker)
      m_faces[face_index].m_marker = f.marker;
    assert( m_faces[face_index].m_vtk_id == f.vtk_id );
  }

  for (const size_t vertex : f.vertices)
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


void Mesh::split_cell(Cell cell, const angem::Plane<double> & plane,
                      const int splitting_face_marker)
{
  if (!m_cells[cell.index()].is_active())
  {
    const auto children = m_cells[cell.index()].immediate_children();
    // make sure this is a cell with hanging nodes
    assert ( children.size() == 1 );
    std::cout << "go down" << std::endl;
    return split_cell(*children[0], plane, splitting_face_marker);
  }
  std::cout << "split " << cell.index() << " (parent "
            << cell.m_parent << ")"<< std::endl;
  assert (cell.is_active());
  // Bookkeeping:
  //  fill polygroup's internal set with the existing vertex coordinates
  // in order to have a map of those to the global vertex indices,
  // which will come in handy when inserting new splitted cells into grid.
  // We can do it because splitting will insert the same vertices plus
  // those that appeared due to plase-face intersection.
  angem::PolyGroup<double> split;
  std::vector<size_t> global_vertex_indices;
  for (const size_t vertex : cell.vertices())
  {
    split.vertices.insert(m_vertices[vertex]);
    global_vertex_indices.push_back(vertex);
  }

  // the actual geometry happens here
  const std::unique_ptr<angem::Polyhedron<double>> polyhedron = cell.polyhedron();
  angem::split(*polyhedron, plane, split,
               constants::marker_below_splitting_plane,
               constants::marker_above_splitting_plane,
               constants::marker_splitting_plane);

  // check we actually split something
  assert( split.vertices.size() > global_vertex_indices.size() );

  // insert new vertices (those that occured due to splitting)
  for (size_t i = global_vertex_indices.size(); i < split.vertices.size(); ++i)
  {
    size_t new_vertex_index = n_vertices();
    // first check that the vertex isn't already present because it's been split
    const size_t split_v_index = m_vertices_from_cell_splitting.find(split.vertices[i]);
    if (split_v_index == m_vertices_from_cell_splitting.size())
    {
      m_vertices_from_cell_splitting.insert( split.vertices[i] );
      m_vertices_from_cell_splitting_indices.push_back(new_vertex_index);
      m_vertices.push_back(split.vertices[i]);
    }
    else
      new_vertex_index = m_vertices_from_cell_splitting_indices[split_v_index];
    global_vertex_indices.push_back(new_vertex_index);
  }

  // map local indices to global
  std::vector<std::vector<size_t>> face_vertex_global_numbering;
  for (size_t i = 0; i < split.polygons.size(); i++)
      face_vertex_global_numbering.push_back(std::move(
          build_global_face_indices_(split.polygons[i], global_vertex_indices)));

  size_t split_face_local_index;  // need this to figure out parent/child faces
  for (size_t i = 0; i < split.polygons.size(); i++)
    if ( split.markers[i] == constants::marker_splitting_plane )
      split_face_local_index = i;

  // make two groups of faces (polyhedra) that will form the new cells
  std::vector<FaceTmpData> tmp_faces(split.polygons.size());
  std::vector<size_t> cell_above_faces, cell_below_faces;
  std::unordered_map<size_t,std::vector<size_t>> cells_to_insert_hanging_nodes;
  for (size_t i = 0; i < split.polygons.size(); i++)
  {
    FaceTmpData & f = tmp_faces[i];
    f.vertices = face_vertex_global_numbering[i];
    f.vtk_id = face_vtk_id_(f.vertices.size());
    std::pair<bool,size_t> face_parent_match =
        determine_face_parent_(face_vertex_global_numbering[i], cell,
                               face_vertex_global_numbering[split_face_local_index]);
    f.parent = face_parent_match.second;

    if (!face_parent_match.first)
    {
      // std::cout << "face neighbor lookup for hanging nodes" << std::endl;
      for (const auto p_cell : face(f.parent).neighbors())
        if (*p_cell != cell)
          cells_to_insert_hanging_nodes[p_cell->index()].push_back(i);
    }

    if ( split.markers[i] == constants::marker_splitting_plane )
      f.marker = splitting_face_marker;
    else
      f.marker = m_faces[f.parent].marker();

    if ( split.markers[i] == constants::marker_below_splitting_plane ||
         split.markers[i] == constants::marker_splitting_plane )
      cell_below_faces.push_back(i);

    if (split.markers[i] == constants::marker_above_splitting_plane ||
        split.markers[i] == constants::marker_splitting_plane)
      cell_above_faces.push_back(i);
  }

  // insert new cells
  const size_t child_cell_index1 = insert_cell_(cell_above_faces, tmp_faces, cell.marker());
  const size_t child_cell_index2 = insert_cell_(cell_below_faces, tmp_faces, cell.marker());
  m_n_split_cells++;

  // handle parent/child cell dependencies
  m_cells[cell.index()].m_children = {child_cell_index1, child_cell_index2};
  m_cells[child_cell_index1].m_parent = cell.index();
  m_cells[child_cell_index2].m_parent = cell.index();

  for (const auto &it : cells_to_insert_hanging_nodes)
  {
    // std::cout << "insert hanging into " << it.first
    //           << "(" << this->cell(it.first).ultimate_parent().index()
    //           << ")"<< std::endl;
    assert( this->cell(it.first).is_active() );
    insert_cell_with_hanging_nodes_(this->cell(it.first), tmp_faces, it.second);
  }
}

active_cell_const_iterator Mesh::begin_active_cells() const
{
  for (auto cell = begin_cells(); cell != end_cells(); ++ cell)
    if (!cell->is_active()) cell++;
    else return active_cell_const_iterator(&*cell);
  return active_cell_const_iterator(nullptr);
}

active_cell_iterator Mesh::begin_active_cells()
{
  for (auto cell = begin_cells(); cell != end_cells(); ++ cell)
    if (!cell->is_active()) cell++;
    else return active_cell_iterator(&*cell);
  return active_cell_iterator(nullptr);
}

std::vector<std::size_t>
Mesh::build_global_face_indices_(const std::vector<size_t> & polygon_local_indices,
                                 const std::vector<size_t> & local_to_global) const
{
  std::vector<std::size_t> global_face_indices(polygon_local_indices.size());
  for (std::size_t i=0; i<polygon_local_indices.size(); ++i)
    global_face_indices[i] = local_to_global[polygon_local_indices[i]];
  return global_face_indices;
}

std::pair<bool,std::size_t> Mesh::determine_face_parent_(const std::vector<size_t> & face_vertices,
                                                         const Cell                & parent_cell,
                                                         const std::vector<size_t> & splitting_face_vertices) const
{
  if ( face_vertices == splitting_face_vertices )
    return {true, constants::invalid_index };  // new face

  // if a face has common vertices with the splitting face then it is a chld face
  std::set<size_t> common_vertices;
  for (const size_t v : face_vertices)
    common_vertices.insert(v);
  for (const size_t v : splitting_face_vertices)
    common_vertices.insert(v);
  if ( common_vertices.size() < face_vertices.size() + splitting_face_vertices.size())
  {
    // there are common vertices + the face normal should match
    for (const Face* face : parent_cell.faces())
    {
      std::set<size_t> common_vertices;
      for (const size_t v : face_vertices)
        common_vertices.insert(v);
      for (const size_t v : face->vertices())
        common_vertices.insert(v);

      if (common_vertices.size() < face_vertices.size() + face->vertices().size())
      { // now compare normals
        const angem::Polygon<double> face_polygon(m_vertices, face_vertices);
        if ( face_polygon.normal().parallel( face->normal(), 1e-4 ) )
          return {false, face->index()};

      }
    }
  }
  else
  {
    // no common vertices: face should match its parent
    std::vector<size_t> sorted = sort_copy_(face_vertices);
    for (const Face* face : parent_cell.faces())
    {
      std::vector<size_t> other_sorted = sort_copy_( face->vertices() );
      if ( sorted == other_sorted )
        return {true, face->index() };
    }
  }

  throw std::runtime_error("Should not be here");
  return {true, constants::invalid_index };  // shut up compiler
}

active_face_const_iterator Mesh::begin_active_faces() const
{
  for (auto face = begin_faces(); face != end_faces(); ++face)
    if (!face->is_active()) face++;
    else return active_face_const_iterator(&*face, m_faces);
  return active_face_const_iterator(nullptr, m_faces);  // end iterator
}

void Mesh::coarsen_cells()
{
  // find all deleted cells
  size_t min_cell_delete_index = std::numeric_limits<size_t>::max();
  for (auto cell = begin_cells(); cell != end_cells(); ++cell)
  {
    if (!cell->m_children.empty())
      cell->m_children.clear();
    else if (cell->parent() != *cell)  // to be deleted
      min_cell_delete_index = std::min(min_cell_delete_index, cell->index());
  }

  //  clear unused vertices
  for ( auto & vertex_cells : m_vertex_cells )
    for (auto it_cell = vertex_cells.begin(); it_cell != vertex_cells.end();)
      if (*it_cell >= min_cell_delete_index)
        vertex_cells.erase(it_cell);
      else ++it_cell;

  // find minimum vertex to erase: new vertices are always at the end
  size_t min_vertex_to_delete = std::numeric_limits<size_t>::max();
  for (std::size_t i=0; i<m_vertex_cells.size(); ++i)
    if (m_vertex_cells[i].empty())
    {
      min_vertex_to_delete = i;
      break;
    }
  m_vertex_cells.erase(m_vertex_cells.begin() + min_vertex_to_delete, m_vertex_cells.end() );

  // clear faces: if face has a deleted vertex then delete it
  // these faces are also consequtive and put into the end
  size_t min_face_to_delete = std::numeric_limits<size_t>::max();
  for (const auto & vertex_faces : m_vertex_faces)
    for (const size_t iface : vertex_faces )
      min_face_to_delete = std::min( min_face_to_delete, iface );
  m_faces.erase( m_faces.begin() + min_face_to_delete, m_faces.end() );
  m_vertex_faces.erase( m_vertex_faces.begin() + min_vertex_to_delete, m_vertex_faces.end() );
  // delete cells
  m_cells.erase( m_cells.begin() + min_cell_delete_index, m_cells.end() );
  // no more split cells
  m_n_split_cells = 0;
  m_n_cells_with_hanging_nodes = 0;
}

bool Mesh::insert_cell_with_hanging_nodes_(Cell & parent,
                                             std::vector<FaceTmpData> big_face_vector,
                                             std::vector<size_t> split_faces)
{
  assert(split_faces.size() == 2);
  const auto & child_face1 = big_face_vector[split_faces[0]];
  const auto & child_face2 = big_face_vector[split_faces[1]];
  assert( child_face1.parent == child_face2.parent );
  const auto c1verts = sort_copy_(child_face1.vertices);
  const auto c2verts = sort_copy_(child_face2.vertices);
  // copy untouched faces
  size_t match_counter = 0;
  for (const auto p_face : parent.faces())
  {
    if (p_face->index() != child_face1.parent)
    {
      const auto sorted_fverts = sort_copy_(p_face->vertices());
      if ( sorted_fverts == c1verts or sorted_fverts == c2verts )
      {
        match_counter++;
        continue;
      }

      FaceTmpData f;
      f.vertices = p_face->vertices();
      f.parent = p_face->parent().index();
      f.marker = p_face->marker();
      f.vtk_id = p_face->vtk_id();
      split_faces.push_back( big_face_vector.size() );
      big_face_vector.push_back(std::move(f));
    }
  }

  if (match_counter == 0)
  {
    const size_t child_cell_index = insert_cell_(split_faces, big_face_vector, parent.marker());
    std::cout << "insert hanging " << child_cell_index << " (" <<
        parent.index() << " ult " << parent.ultimate_parent().index()
              << ")"<< std::endl;
    m_cells[child_cell_index].m_parent = parent.index();
    m_cells[parent.index()].m_children = {child_cell_index};
    m_n_cells_with_hanging_nodes++;
    return true;
  }
  else assert ( match_counter == 2 );
  return false;
}

// const Face & Mesh::ultimate_parent(const Face & face) const
// {
//   const Face * par = &m_faces[ face.parent() ];
//   while (par->parent() != par->index())
//     par = &m_faces[ par->parent() ];
//   return *par;
// }

// Face & Mesh::ultimate_parent(const Face & face)
// {
//   return const_cast<Face&>(ultimate_parent(face));
// }

}  // end namespace mesh
