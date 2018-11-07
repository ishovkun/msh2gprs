#include <Mesh.hpp>

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
                  const int type)
{
  if (type == -1)
  {
    const int n_verts = poly.get_points().size();
    // if (n_verts == 4)
      // shape_ids.push_back();
  }

  std::vector<std::size_t> indices;
  const std::vector<Point> & points = poly.get_points();
  for (const auto & p : points)
  {
    const std::size_t ind = vertices.insert(p);
    indices.push_back(ind);
  }

  cells.push_back(indices);

  for (const auto & edge : poly.get_faces())
  {

  }


}

}
