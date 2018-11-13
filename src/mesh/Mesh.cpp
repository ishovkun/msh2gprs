#include <Mesh.hpp>

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
      (iter->second).push_back(new_element_index);
    else
      map_faces.insert({ {hash, {new_element_index}} });
  }
}


const std::vector<std::size_t> & Mesh::get_neighbors( const Face & face ) const
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

  return iter->second;
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
  // find  global indices of polygon vertices
  std::vector<std::size_t> indices;
  const std::vector<Point> & points = poly.get_points();
  for (const auto & p : points)
  {
    const std::size_t ind = vertices.find(p);
    indices.push_back(ind);
  }

  // get faces with global indexing of vertices
  std::vector<std::vector<std::size_t>> faces;
  for (const auto & face : poly.get_faces())
  {
    std::vector<std::size_t> face_glob;
    for (const auto ivert : face)
      face_glob.push_back(indices[ivert]);
    faces.push_back(face_glob);
  }
  return std::move(faces);

}


void Mesh::insert_physical_face(const Polygon & poly,
                                const int       marker)
{
  const auto & points = poly.get_points();
  std::vector<std::size_t> face(points.size());
  for (int i=0; i<points.size(); ++i)
    face[i] = vertices.find(points[i]);

  const auto hash = hash_value(face);
  auto it = map_physical_faces.find(hash);
  if (it == map_physical_faces.end())
    map_physical_faces.insert({hash, marker});
  else
    it->second = marker;
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


}
