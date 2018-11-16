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
    {
      if (jcell != icell)
      {
        neighbors.push_back(jcell);
      }
    }
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


// void Mesh::split(face_iterator & face)
// {
//   auto new_points = face.vertices();
//   const std::size_t old_n_vertices = vertices.size();
//   std::vector<std::size_t> indices;
//   indices.reserve(new_points.size());
//   std::size_t index = old_n_vertices;
//   for (auto vertex : vertices)
//   {
//     vertices.points.push_back(vertex);
//     indices.push_back(index);
//     index++;
//   }

//   Face new_face;
//   new_face.index = n_faces();
//   new_face.marker = face.marker();
//   const auto hash = hash_value(indices);

//   std::cout << "not implemented" << std::endl;
//   exit(0);
//   // TODO: update element vertices
//   // duplication of vertices depends on how many fractures
//   // intersect in a vertex
//   // auto neighbors = face.neighbors();
//   // if (neighbors.size() > 1)
//   // {
//   //   for (auto & neighbor : neighbors)
//   //   {

//   //   }
//   }
// }


void Mesh::mark_for_split(const face_iterator & face)
{
  marked_for_split.push_back(face.hash());
}


void Mesh::split_faces()
{
  /* Algorithm:
  * 1. create SurfaceMesh from marked faces
  * 2. find internal vertices (not on boundary)
  * 3. split each face that has internal vertices
  */
  SurfaceMesh<double> mesh_faces(1e-6);
  for (const auto & hash : marked_for_split)
  {
    face_iterator face = create_face_iterator(map_faces.find(hash));
    Polygon poly(face.vertices());
    mesh_faces.insert(poly);
  }

  std::unordered_set<std::size_t> internal_points =
      find_internal_vertices(mesh_faces);

}

}
