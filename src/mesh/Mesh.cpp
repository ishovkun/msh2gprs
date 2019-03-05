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
      face_data.old_index = face_data.index;
      map_faces.insert({ {hash, face_data} });
    }
  }
}


void Mesh::insert_cell(const std::vector<std::size_t> & ivertices,
                       const int                        vtk_id,
                       const int                        marker)
{
  const std::vector<std::vector<std::size_t>> poly_faces =
      angem::PolyhedronFactory::get_global_faces<double>(ivertices, vtk_id);

  const std::size_t new_element_index = cells.size();

  cells.push_back(ivertices);
  shape_ids.push_back(vtk_id);
  cell_markers.push_back(marker);

  for (const auto & face : poly_faces)
  {
    const auto hash = hash_value(face);
    auto iter = map_faces.find(hash);
    if (iter != map_faces.end())
      (iter->second).neighbors.push_back(new_element_index);
    else
    {
      Face face_data;
      face_data.neighbors.push_back(new_element_index);
      face_data.index = map_faces.size();
      face_data.old_index = face_data.index;
      face_data.ordered_indices = face;
      switch (face.size())
      {
        case 3:
          face_data.vtk_id = 5;  //  triangle
          break;
        case 4:
          face_data.vtk_id = 9;  //  vtk_quad
          break;
        default:
          face_data.vtk_id = 7;  //  vtk_polygon
          break;
      }
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
  const std::vector<std::vector<std::size_t>> cell_faces =
      angem::PolyhedronFactory::get_global_faces<double>(cells[icell],
                                                         shape_ids[icell]);
  for (const auto & face : cell_faces)
  {
    const std::vector<std::size_t> & face_neighbors = get_neighbors(face);
    for (const std::size_t jcell : face_neighbors)
      if (jcell != icell)
        neighbors.push_back(jcell);
  }
  return neighbors;
}


std::vector<std::vector<std::size_t>> Mesh::get_faces(const Polyhedron & poly) const
{
  return get_face_indices(poly, vertices);
}


void Mesh::insert(const Polygon & poly,
                  const int       marker)
{
  const auto & points = poly.get_points();
  std::vector<std::size_t> face(points.size());
  for (int i=0; i<points.size(); ++i)
    face[i] = vertices.find(points[i]);

  // const vtk_id = poly.vtk_id();
  const int vtk_id = -1;
  insert_face(face, vtk_id, marker);
}


void Mesh::insert_face(const std::vector<std::size_t> & ivertices,
                       const int                        vtk_id,
                       const int                        marker)
{
  const auto hash = hash_value(ivertices);
  auto it = map_faces.find(hash);
  if (it == map_faces.end())
  {
    Face face_data;
    face_data.marker = marker;
    face_data.vtk_id = vtk_id;
    face_data.index = map_faces.size();
    face_data.old_index = face_data.index;
    face_data.ordered_indices = ivertices;
    map_faces.insert({hash, face_data});
  }
  else
  {
    it->second.marker = marker;
    it->second.vtk_id = vtk_id;
  }
}


Point Mesh::get_center(const std::size_t icell) const
{
  return get_element_center(vertices, cells[icell]);
}


std::unique_ptr<Polyhedron> Mesh::get_polyhedron(const std::size_t icell) const
{
  return angem::PolyhedronFactory::create<double>(vertices.points,
                                                  cells[icell],
                                                  shape_ids[icell]);
}


void Mesh::split_vertex(const std::size_t                              ivertex,
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

  std::vector<std::vector<std::size_t>> groups =
      group_cells_based_on_split_faces(affected_cells, vertex_faces);
  const int n_groups = groups.size();

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
}


void Mesh::mark_for_split(const face_iterator & face)
{
  marked_for_split.push_back(face.hash());
}


SurfaceMesh<double> Mesh::split_faces()
{
  // std::cout << "n_nodes() = " << n_vertices() << std::endl;
  // std::cout << "n_faces = " << n_faces() << std::endl;
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

  const std::size_t n_faces_old = this->n_faces();
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
        auto it = create_face_iterator
            (map_faces.find(map_2d_3d.find(neighbor)->second));
        vertex_faces.push_back(std::move(it));
      }

      const auto edge_vertices = edge.vertices();
      std::cout << "splitting vertex = " << edge_vertices.first << std::endl;
      split_vertex(vertices.find(edge_vertices.first), vertex_faces,
                   map_old_new_cells, new_cells);
      std::cout << "splitting vertex = " << edge_vertices.second << std::endl;
      split_vertex(vertices.find(edge_vertices.second), vertex_faces,
                   map_old_new_cells, new_cells);
    }
  }

  std::cout << "modifying face map" << std::endl;

  // modify face map
  std::unordered_set<std::size_t> old_ind_touched;
  std::unordered_set<hash_type> faces_to_delete;
  // store it since new faces are added first and then old faces are deleted
  // which may cause wrong indexing
  std::size_t num_faces = map_faces.size();
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

    assert(old_poly_faces.size() == new_poly_faces.size());

    for (std::size_t i=0; i<old_poly_faces.size(); ++i)
    {
      const hash_type old_hash = hash_value(old_poly_faces[i]);
      const hash_type new_hash = hash_value(new_poly_faces[i]);

      if (new_hash != old_hash) // face changed
      {
        auto it_face = map_faces.find(new_hash);
        auto it_face_old = map_faces.find(old_hash);
        if (it_face == map_faces.end())
        {
          Face new_face;
          new_face.neighbors.push_back(icell);
          new_face.marker = it_face_old->second.marker;
          new_face.ordered_indices = new_poly_faces[i];

          if (old_ind_touched.insert(it_face_old->second.index).second) // if not touched
            new_face.index = it_face_old->second.index;
          else
          {
            new_face.index = num_faces;
            num_faces++;
          }

          new_face.old_index = it_face_old->second.old_index;
          new_face.vtk_id = it_face_old->second.vtk_id;
          map_faces.insert({new_hash, new_face});
        }
        else
          it_face->second.neighbors.push_back(icell);

        // mark old faces for delete
        faces_to_delete.insert(old_hash);
      }  // end if face changed
    }    // end face loop

    // replace old cell with new cell
    cells[icell] = new_cells[new_icell];
  }

  // remove marked old faces from map
  for (auto & hash : faces_to_delete)
  {
    auto face_it = map_faces.find(hash);
    if (face_it != map_faces.end())
      map_faces.erase(face_it);
  }

  // clear marked elements vector
  marked_for_split.clear();
  return mesh_faces;
}


std::vector<std::vector<std::size_t>>
Mesh::group_cells_based_on_split_faces(const std::unordered_set<std::size_t> & affected_cells,
                                       const std::vector<face_iterator>      & vertex_faces) const
{
  // group affected elements
  // two elements are in the same group if they are neighbors and
  // the neighboring face is not in vertex_faces array
  const int n_groups = vertex_faces.size();
  std::unordered_map<std::size_t, int> map_cell_group;
  int igroup = 0;
  int new_group = 0;
  // just is purely for faster checking: cells that are already processed
  std::unordered_set<std::size_t> processed_cells;
  for (const std::size_t icell : affected_cells)
  {
    // std::cout << "icell = " << icell << std::endl;

    auto group_it = map_cell_group.find(icell);
    if (group_it != map_cell_group.end())
    {
      igroup = group_it->second;
    }
    else
    {
      igroup = new_group;
      // groups[igroup].push_back(icell);
      map_cell_group.insert({icell, igroup});
      new_group++;
    }

    processed_cells.insert(icell);

    for (const std::size_t jcell : get_neighbors(icell))
      // if (processed_cells.find(jcell) == processed_cells.end() and
      //     affected_cells.find(jcell) != affected_cells.end())
      if (affected_cells.find(jcell) != affected_cells.end())
      {
        auto pair_cells = std::minmax(icell, jcell);
        std::vector<std::size_t> ordered_neighbors =
            {pair_cells.first, pair_cells.second};

        // find out if i and j neighbor by a marked face
        bool neighbor_by_marked_face = false;
        for (const auto & face : vertex_faces)
        {
          const auto it = map_faces.find(face.hash());

          if (it->second.neighbors == ordered_neighbors)
          {
            neighbor_by_marked_face = true;
            break;
          }
        }

        if (!neighbor_by_marked_face)
        {
          // std::cout << " same" << std::endl;
          auto group_it = map_cell_group.find(jcell);
          if (group_it == map_cell_group.end())
            map_cell_group.insert({jcell, igroup});
          else
          {
            if (group_it->second < igroup)
            {
              map_cell_group[icell] = group_it->second;
              igroup = group_it->second;
              new_group--;
            }
            else
              map_cell_group[jcell] = igroup;
          }
        }
        // else
        //   std::cout << " other" << std::endl;
      }
    // std::cout << "igroup = " << igroup << std::endl << std::endl;
  }
  // std::cout << "groups" << std::endl;
  // for (auto & it : map_cell_group)
  //   std::cout << it.first << "\t" << it.second << std::endl;

  // abort();

  std::vector<std::vector<std::size_t>> groups(n_groups);
  for (auto it : map_cell_group)
    groups[it.second].push_back(it.first);

  return groups;
}


std::vector<face_iterator> Mesh::get_ordered_faces()
{
  std::vector<face_iterator> ordered_faces(this->n_faces(), this->begin_faces());
  for (auto face=begin_faces(); face!=end_faces(); ++face)
    ordered_faces[face.index()] = face;
  return ordered_faces;
}

}
