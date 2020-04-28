#include "FaceSplitter.hpp"

namespace mesh {

FaceSplitter::FaceSplitter(const Mesh & grid)
    : _grid(grid), _vertex_coord(_grid.vertices())
{
  _cell_vertices.resize(_grid.n_cells());
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
    _cell_vertices[cell->index()] = cell->vertices();
}

void FaceSplitter::mark_for_split(const size_t face_index)
{
  if (!_grid.face(face_index).is_active())
    throw std::invalid_argument("cannot split an inactive face");

  _marked_for_split.push_back(face_index);
}

SurfaceMesh<double> FaceSplitter::split_faces()
{
  // create surfacemesh and map vertices
  SurfaceMesh<double> mesh_faces(1e-6);
  // map 2d-element -> 3d face hash
  std::unordered_map<std::size_t, std::size_t> map_face_surface_element;
  // map surfacemesh vertex -> 3d mesh vertex
  std::unordered_map<std::size_t, std::size_t> map_frac_vertex_vertex;

  /* Algorithm:
  * create SurfaceMesh from marked faces in order to identify
  * vertices to split those whose edge have >1 neighbors)
  * cross-match vertices in 3d Mesh and Surface mesh. */
  create_fracture_face_grid_(mesh_faces, map_face_surface_element, map_frac_vertex_vertex);

  // Identify the vertices that require splitting
  std::unordered_map<size_t, std::vector<size_t>> vertices_to_split;
  find_vertices_to_split_(mesh_faces, map_face_surface_element,
                          map_frac_vertex_vertex, vertices_to_split);

  // split 'em
  for (const auto & it : vertices_to_split)
  {
    const size_t vertex = it.first;
    const auto & faces = it.second;
    std::cout << "splitting vertex = " << vertex << std::endl;
    split_vertex_(vertex, faces);
  }

  return mesh_faces;
}

void FaceSplitter::split_vertex_(const std::size_t               vertex_index,
                                 const std::vector<std::size_t> &splitted_face_indices)
{
  std::vector<std::vector<std::size_t>> groups =
      group_cells_based_on_split_faces_(_grid.m_vertex_cells[vertex_index], splitted_face_indices);

  // create new vertices
  std::vector<std::size_t> new_vertex_indices(groups.size());
  const angem::Point<3,double> vertex_coord = _vertex_coord[vertex_index];
  for (std::size_t group = 0; group < groups.size(); group++)
  {
    if (group == 0)  // group 0 retains old vertex
      new_vertex_indices[group] = vertex_index;
    else  // add new vertices
    {
      const std::size_t new_vertex_index = _vertex_coord.size();
      new_vertex_indices[group] = new_vertex_index;
      _vertex_coord.push_back(vertex_coord);
    }
  }

  // modify cell vertices: replace vertex indices with the new vertices
  // start from 1 since group 0 retains the old index
  for (std::size_t group = 1; group < groups.size(); group++)
  {
    const std::vector<size_t> & cell_group = groups[group];
    for (const std::size_t cell_index : cell_group)
    {
      std::vector<size_t> & cell_vertices = _cell_vertices[cell_index];
      for (size_t & cell_vertex_index : cell_vertices)
        if (cell_vertex_index == vertex_index)
        {
          cell_vertex_index = new_vertex_indices[group];
        }
    }
  }
}

std::vector<std::vector<std::size_t>>
FaceSplitter::
group_cells_based_on_split_faces_(const std::vector<size_t> & affected_cells,
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
    for (const Cell* jcell : _grid.cell(icell).neighbors())
    {
      if (std::find(affected_cells.begin(), affected_cells.end(),
                    jcell->index()) != affected_cells.end())
      {
        // take index explicitly since minmax takes a reference
        const size_t jind = jcell->index();
        // what face neighbors should be
        const auto pair_cells = std::minmax(icell, jind);

        // find out if cell i and cell j neighbor by a marked face
        bool neighbor_by_marked_face = false;
        for (const size_t iface : split_faces)
        {
          const Face f = _grid.face(iface);
          const auto f_neighbors = f.neighbors();

          assert( f_neighbors.size() == 2 );
          auto pair_cells2 = std::minmax(f_neighbors[0]->index(),
                                         f_neighbors[1]->index());

          if (pair_cells == pair_cells2)
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

void FaceSplitter::
find_vertices_to_split_(const SurfaceMesh<double> & mesh_faces,
                        const std::unordered_map<std::size_t, std::size_t> & map_face_surface_element,
                        std::unordered_map<std::size_t,std::size_t> & map_frac_vertex_vertex,
                        std::unordered_map<size_t, std::vector<size_t>> & vertices_to_split) const
{
  for (auto edge = mesh_faces.begin_edges(); edge !=mesh_faces.end_edges(); ++edge)
  {
    const std::vector<size_t> & edge_neighbors = edge.neighbors();
    std::vector<size_t> grid_face_indices;
    grid_face_indices.reserve(edge_neighbors.size());
    // if edge is a neighbor of a split face, then consider it
    for (size_t ielement: edge_neighbors)
    {
      const auto it = map_face_surface_element.find(ielement);
      if (it == map_face_surface_element.end())
        throw std::runtime_error("error during splitting faces");

      grid_face_indices.push_back(it->second);
    }

    // add edge vertices to the map vertex->vector_of_affected_faces
    if (edge_neighbors.size() > 1)  // internal edge vertex
    {
      const auto edge_vertices = edge.vertex_indices();
      const size_t v1 = map_frac_vertex_vertex[edge_vertices.first];
      const size_t v2 = map_frac_vertex_vertex[edge_vertices.second];
      auto it1 = vertices_to_split.find(v1);
      if (it1 == vertices_to_split.end())
        vertices_to_split.insert({v1, grid_face_indices});
      else
        for (const size_t face : grid_face_indices)
          if ( std::find(it1->second.begin(), it1->second.end(), face) == it1->second.end())
            it1->second.push_back(face);

      auto it2 = vertices_to_split.find(v2);
      if (it2 == vertices_to_split.end())
        vertices_to_split.insert({v2, grid_face_indices});
      else
        for (const size_t face : grid_face_indices)
          if ( std::find(it2->second.begin(), it2->second.end(), face) ==
               it2->second.end())
            it2->second.push_back(face);
    }
  }
}

void FaceSplitter::create_fracture_face_grid_(SurfaceMesh<double> & mesh_faces,
                                  std::unordered_map<size_t,size_t> & map_face_surface_element,
                                  std::unordered_map<size_t,size_t> & map_frac_vertex_vertex) const
{
  for (const size_t face_index: _marked_for_split)
  {
    const Face & f = _grid.face(face_index);
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
}

}  // end namespace mesh
