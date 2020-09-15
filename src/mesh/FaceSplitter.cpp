#include "FaceSplitter.hpp"
#include "FractureTip.hpp"
#include "Groups.hpp"


namespace mesh {

FaceSplitter::FaceSplitter(const Mesh & grid)
    : _grid(grid), _vertex_coord(_grid.vertices())
{
  _cell_vertices.resize(_grid.n_cells_total());
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
  /* Algorithm:
  * create SurfaceMesh from marked faces in order to identify
  * vertices to split those whose edge have >1 neighbors)
  * cross-match vertices in 3d Mesh and Surface mesh. */
  SurfaceMesh<double> mesh_faces(1e-6);
  create_fracture_face_grid_(mesh_faces, _surface_to_face, _surface_vertex_to_global);

  // Identify the vertices that require splitting
  find_vertices_to_split_(mesh_faces);

  // split 'em
  for (const auto & it : _vertices_to_split)
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
  std::vector<size_t> affected_active_cells;
  std::back_insert_iterator< std::vector<size_t> > back_it (affected_active_cells);
  std::copy_if( _grid.m_vertex_cells[vertex_index].begin(), _grid.m_vertex_cells[vertex_index].end(),
                back_it, [this](const size_t icell){ return _grid.cell(icell).is_active();});
  std::vector<std::vector<std::size_t>> groups =
      group_cells_based_on_split_faces_2(affected_active_cells, splitted_face_indices);
  // std::vector<std::vector<std::size_t>> groups =
  //     group_cells_based_on_split_faces_2(_grid.m_vertex_cells[vertex_index], splitted_face_indices);
  // assert( groups.size() > 1 );
  if ( groups.size() == 1 )
    return;

  // create new vertices
  std::vector<std::size_t> new_vertex_indices(groups.size());
  auto & child_vertices = _parent_to_child_vertices[vertex_index];
  const angem::Point<3,double> vertex_coord = _vertex_coord[vertex_index];
  for (std::size_t group = 0; group < groups.size(); group++)
  {
    if (group == 0)  // group 0 retains old vertex
      new_vertex_indices[group] = vertex_index;
    else  // add new vertices
    {
      const std::size_t new_vertex_index = _vertex_coord.size();
      // std::cout << "adding vertex " << new_vertex_index << std::endl;
      new_vertex_indices[group] = new_vertex_index;
      _vertex_coord.push_back(vertex_coord);
      child_vertices.push_back(new_vertex_index);
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
group_cells_based_on_split_faces_2(const std::vector<size_t> & affected_cells,
                                   const std::vector<size_t> & split_faces) const
{
  const size_t n_groups = std::max(split_faces.size(), size_t(2));
  utils::Groups groups(affected_cells);
  int new_group = 0;
  for (const std::size_t icell : affected_cells)
    if (_grid.cell(icell).is_active())
    {
      for (auto face : _grid.cell(icell).faces() )
        if (std::find(split_faces.begin(), split_faces.end(), face->index()) == split_faces.end())
        {
          const auto face_neighbors = face->neighbors();

          if (face_neighbors.size() == 1)
            continue;

          size_t jcell = face_neighbors[0]->index();
          if (jcell == icell)
            jcell = face_neighbors[1]->index();

          if (std::find(affected_cells.begin(), affected_cells.end(), jcell) != affected_cells.end())
          {
            assert( _grid.cell(jcell).is_active() );
            groups.merge(icell, jcell);
          }
        }
    }

  return groups.get();
}

void FaceSplitter::
find_vertices_to_split_(const SurfaceMesh<double> & mesh_faces)
{
  mesh::FractureTip tip(mesh_faces, _grid, _surface_vertex_to_global);
  for (auto edge = mesh_faces.begin_edges(); edge !=mesh_faces.end_edges(); ++edge)
  {
    const std::vector<size_t> & edge_neighbors = edge.neighbors();
    std::vector<size_t> grid_face_indices;
    grid_face_indices.reserve(edge_neighbors.size());

    // if edge is a neighbor of a split face, then consider it
    for (size_t ielement: edge_neighbors)
    {
      const auto it = _surface_to_face.find(ielement);
      if (it == _surface_to_face.end())
        throw std::runtime_error("error during splitting faces");

      grid_face_indices.push_back(it->second);
    }

    // add edge vertices to the map vertex->vector_of_affected_faces
    if (edge_neighbors.size() > 1)  // internal edge vertex
    {
      const auto edge_vertices = edge.vertex_indices();
      const size_t v1 = _surface_vertex_to_global[edge_vertices.first];
      const size_t v2 = _surface_vertex_to_global[edge_vertices.second];
      if (!tip.contains(v1))
      {
        auto it1 = _vertices_to_split.find(v1);
        if (it1 == _vertices_to_split.end())
          _vertices_to_split.insert({v1, grid_face_indices});
        else
          for (const size_t face : grid_face_indices)
            if (std::find(it1->second.begin(), it1->second.end(), face) == it1->second.end())
              it1->second.push_back(face);
      }

      if (!tip.contains(v2))
      {
        auto it2 = _vertices_to_split.find(v2);
        if (it2 == _vertices_to_split.end())
          _vertices_to_split.insert({v2, grid_face_indices});
        else
          for (const size_t face : grid_face_indices)
            if ( std::find(it2->second.begin(), it2->second.end(), face) == it2->second.end())
              it2->second.push_back(face);
      }
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
