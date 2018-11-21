#include <Mesh.hpp>
#include <SurfaceMesh.hpp>

#include <angem/PolyhedronFactory.hpp>

namespace mesh
{


Mesh::Mesh()
{}


void Mesh::insert(const Polyhedron & poly,
                  const int          marker)
{
  std::vector<std::size_t> indices;
  const std::vector<Point> & points = poly.get_points();
  for (const auto & p : points)
  {
    const std::size_t ind = vertices.insert(p);
    indices.push_back(ind);
  }

  const std::size_t new_element_index = cells.size();
  cells.push_back(indices);
  shape_ids.push_back(poly.id());
  cell_markers.push_back(marker);

  for (const auto & face : poly.get_faces())
  {
    // get global indices of face vertices
    std::vector<std::size_t> face_glob;
    for (const auto & ivert : face)
      face_glob.push_back(indices[ivert]);

    const auto hash = hash_value(face_glob);

    auto iter = map_faces.find(hash);
    if (iter != map_faces.end())
    {
      (iter->second).neighbors.push_back(new_element_index);
    }
    else
    {
      Face face_data;
      face_data.neighbors.push_back(new_element_index);
      face_data.index = map_faces.size();
      map_faces.insert({ {hash, face_data} });
    }
  }
}


const std::vector<std::size_t> & Mesh::get_neighbors( const FaceiVertices & face ) const
{
  const auto hash = hash_value(face);
  const auto iter = map_faces.find(hash);

  if (iter == map_faces.end())
  {
    for (auto & v : face)
      std::cout << v << " ";
    std::cout << std::endl;
    throw std::out_of_range("face does not exist");
  }

  return iter->second.neighbors;
}


std::vector<std::size_t>
Mesh::get_neighbors( const std::size_t icell ) const
{
  if (icell >= cells.size())
    throw std::out_of_range("wrong cell index: " + std::to_string(icell));

  std::vector<std::size_t> neighbors;
  const Polyhedron poly =
      angem::PolyhedronFactory::create<double>(vertices.points, cells[icell],
                                               shape_ids[icell]);
  for (const auto & face : get_faces(poly))
  {
    const std::vector<std::size_t> & face_neighbors = get_neighbors(face);
    for (const std::size_t jcell : face_neighbors)
      if (jcell != icell)
        neighbors.push_back(jcell);
  }
  return std::move(neighbors);
}


std::vector<std::vector<std::size_t>> Mesh::get_faces(const Polyhedron & poly) const
{
  return std::move(get_face_indices(poly, vertices));
}


void Mesh::insert(const Polygon & poly,
                  const int       marker)
{
  const auto & points = poly.get_points();
  std::vector<std::size_t> face(points.size());
  for (int i=0; i<points.size(); ++i)
    face[i] = vertices.find(points[i]);

  const auto hash = hash_value(face);
  auto it = map_faces.find(hash);
  if (it == map_faces.end())
  {
    Face face_data;
    face_data.marker = marker;
    face_data.index = map_faces.size();
    map_faces.insert({hash, face_data});
  }
  else
    it->second.marker = marker;
}


Point Mesh::get_center(const std::size_t icell) const
{
  return std::move(get_element_center(vertices, cells[icell]));
}


Polyhedron Mesh::get_polyhedron(const std::size_t icell) const
{
  return angem::PolyhedronFactory::create<double>(vertices.points,
                                                  cells[icell],
                                                  shape_ids[icell]);
}


void Mesh::split_vertex(const std::size_t                            ivertex,
                        const std::vector<face_iterator>             & vertex_faces,
                        std::unordered_map<std::size_t, std::size_t> & map_old_new_cells,
                        std::vector<std::vector<std::size_t>>        & new_cells)
{
  // find affected elements
  std::unordered_set<std::size_t> affected_cells;
  for (auto & face : vertex_faces)
    for (const std::size_t icell : face.neighbors())
    {
      affected_cells.insert(icell);
      // include elements that don't neighbor split faces (only by vertex)
      for (const std::size_t jcell : get_neighbors(icell))
      {
        auto cell_j = create_cell_iterator(jcell);
        if (cell_j.has_vertex(ivertex))
          affected_cells.insert(cell_j.index());
      }
    }

  // group affected elements
  // two elements are in the same group if they are neighbors and
  // the neighboring face is not in vertex_faces array
  const int n_groups = vertex_faces.size();
  std::vector<std::vector<std::size_t>> groups(n_groups);
  int igroup = 0;
  // just is purely for faster checking
  std::unordered_set<std::size_t> grouped_cells;
  for (const std::size_t icell : affected_cells)
  {
    if (grouped_cells.find(icell) == grouped_cells.end())
    {
      grouped_cells.insert(icell);
      groups[igroup].push_back(icell);
    }
    else
      continue;

    for (const std::size_t jcell : get_neighbors(icell))
    {
      std::vector<std::size_t> pair_cells;
      if (icell > jcell)
      {
        pair_cells.push_back(jcell);
        pair_cells.push_back(icell);
      }
      else  // they are always not equal
      {
        pair_cells.push_back(icell);
        pair_cells.push_back(jcell);
      }

      // find out if i and j neighbor by a marked face
      bool neighbor_by_marked_face = false;
      for (const auto & face : vertex_faces)
        if (map_faces.find(face.hash())->second.neighbors == pair_cells)
        {
          neighbor_by_marked_face = true;
          break;
        }

      if (!neighbor_by_marked_face)
        groups[igroup].push_back(jcell);
    }

    igroup++;
  }

  // create new vertices
  std::vector<std::size_t> new_ivertices(n_groups);
  const angem::Point<3,double> vertex = vertices[ivertex];
  for (int i=0; i < n_groups; ++i)
  {
    if (i == 0)  // group 0 retains old vertex
      new_ivertices[i] = ivertex;
    else
    {
      const std::size_t new_ivertex = vertices.size();
      new_ivertices[i] = new_ivertex;
      vertices.points.push_back(vertex);
    }
  }

  // std::cout << "modify cell elements" << std::endl;
  // modify new cell vertices
  for (int igroup = 0; igroup < groups.size(); ++igroup)
  {
    const auto & group = groups[igroup];
    for (const std::size_t icell : group)
    {
      // modify cell vertices
      std::vector<std::size_t> * p_new_cell;
      auto it = map_old_new_cells.find(icell);
      if (it != map_old_new_cells.end())  // cell modified before
        p_new_cell = &(new_cells[it->second]);
      else
      {
        map_old_new_cells.insert({icell, new_cells.size()});
        new_cells.push_back(cells[icell]);
        p_new_cell = &(new_cells.back());
      }

      auto & new_cell = *p_new_cell;
      for (std::size_t & jvertex : new_cell)
        if (jvertex == ivertex)
          jvertex = new_ivertices[igroup];
    }
  }
  // std::cout << "end vertex split" << std::endl;

}


void Mesh::mark_for_split(const face_iterator & face)
{
  marked_for_split.push_back(face.hash());
}


void Mesh::split_faces()
{
  /* Algorithm:
  * 1. create SurfaceMesh from marked faces
  * 2. find internal vertices (not on boundary)
  * 3. split each internal vertex
  * 4. modify neighbors map */

  SurfaceMesh<double> mesh_faces(1e-6);
  // map 2d-element : 3d face hash
  std::unordered_map<std::size_t, hash_type> map_2d_3d;
  for (const auto & hash : marked_for_split)
  {
    face_iterator face = create_face_iterator(map_faces.find(hash));
    Polygon poly(face.vertices());
    const std::size_t ielement = mesh_faces.insert(poly);
    map_2d_3d.insert({ielement, hash});
  }

  std::unordered_map<std::size_t, std::size_t> map_old_new_cells;
  std::vector<std::vector<std::size_t>> new_cells;
  std::unordered_map<hash_type, Face> map_new_faces;
  for (auto edge = mesh_faces.begin_edges(); edge !=mesh_faces.end_edges(); ++edge)
  {
    const std::vector<std::size_t> & edge_neighbors = edge.neighbors();
    if (edge_neighbors.size() > 1)  // internal edge vertex
    {
      // find face elements of 3d mesh corresponding to 2d mesh elements
      std::vector<face_iterator> vertex_faces;
      for (const auto & neighbor : edge_neighbors)
      {
        auto it = create_face_iterator(map_faces.find(map_2d_3d.find(neighbor)->second));
        vertex_faces.push_back(std::move(it));
      }

      const auto edge_vertices = edge.vertices();
      split_vertex(vertices.find(edge_vertices.first), vertex_faces,
                   map_old_new_cells, new_cells);
      split_vertex(vertices.find(edge_vertices.second), vertex_faces,
                   map_old_new_cells, new_cells);
    }
  }

  // std::cout << "modifying face map" << std::endl;

  // modify face map
  for (const auto it_cell : map_old_new_cells)
  {
    const std::size_t icell = it_cell.first;
    const std::size_t new_icell = it_cell.second;  // index in new array

    const std::vector<std::vector<std::size_t>> old_poly_faces =
        angem::PolyhedronFactory::get_global_faces<double>(cells[icell],
                                                           shape_ids[icell]);
    const std::vector<std::vector<std::size_t>> new_poly_faces =
        angem::PolyhedronFactory::get_global_faces<double>(new_cells[new_icell],
                                                           shape_ids[icell]);
    // std::cout << "got some faces" << std::endl;

    assert(old_poly_faces.size() == new_poly_faces.size());

    for (int i=0; i<old_poly_faces.size(); ++i)
    {
      const hash_type old_hash = hash_value(old_poly_faces[i]);
      const hash_type new_hash = hash_value(new_poly_faces[i]);
      if (new_hash == old_hash)
      {
        // face did not change:
        continue;
      }
      else  // face changed
      {
        auto it_face = map_faces.find(new_hash);
        auto it_face_old = map_faces.find(old_hash);
        if (it_face == map_faces.end())
        {
          Face new_face;
          new_face.neighbors.push_back(icell);
          new_face.marker = it_face_old->second.marker;
          map_faces.insert({new_hash, new_face});
        }
        else
        {
          it_face->second.neighbors.push_back(icell);
        }

        // remove old connection
        auto & old_neighbors = it_face_old->second.neighbors;
        for (int j=0; j<old_neighbors.size(); ++j)
          if (old_neighbors[j] == icell)
            old_neighbors.erase(old_neighbors.begin() + j);

        if (old_neighbors.size() == 0)
          map_faces.erase(it_face_old);
      }  // end if face changed
    }    // end face loop

    // replace old cell with new cell
    cells[icell] = new_cells[new_icell];
  }

  // compute new face indices
  std::size_t face_index = 0;
  for (auto & face : map_faces)
  {
    face.second.index = face_index;
    face_index++;
  }

  // clear marked elements vector
  marked_for_split.clear();
}


}
