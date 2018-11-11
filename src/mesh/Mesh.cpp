#include <Mesh.hpp>

#include <angem/PolyhedronFactory.hpp>

namespace mesh
{


std::size_t estimate_max_vertices(const int n_polygon_vertices)
{
  const std::size_t max_vertices =
      static_cast<std::size_t>(
          pow(std::pow(2, 256) - 1, 1.0 / n_polygon_vertices));
  return max_vertices;

}


Mesh::Mesh()
    :
    max_vertices(estimate_max_vertices(max_polygon_size_vertices))
{}


uint256_t Mesh::hash_value(const Face & face) const
{
  Face face_sorted = face;
  std::sort(face_sorted.begin(), face_sorted.end());

  uint256_t mult = 1;
  uint256_t result = 0;
  for (int i=0; i<max_polygon_size_vertices; ++i)
  {
    // hashing starts with 1 since
    // (1,2,3) and (0,1,2,3) are different elementes
    result += 1 + mult * face_sorted[i];
    mult *= max_vertices;
  }
  return result;
}


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

  for (const auto & face : poly.get_faces())
  {
    const auto hash = hash_value(face);
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
    throw std::out_of_range("edge does not exist");

  return iter->second;
}


const std::vector<std::size_t> &
Mesh::get_neighbors( const std::size_t icell ) const
{
  std::vector<std::size_t> v_neighbors;
  const Polyhedron poly =
      angem::PolyhedronFactory::create<double>(vertices.points, cells[icell]);
  for (const auto face : poly.get_faces())
  {
    const std::vector<std::size_t> & face_neighbors = get_neighbors(face);
    for (const std::size_t jcell : face_neighbors)
      if (jcell != icell)
        v_neighbors.push_back(jcell);
  }
  return std::move(v_neighbors);
}


}
