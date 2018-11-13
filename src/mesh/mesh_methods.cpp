#include <mesh_methods.hpp>
#include <angem/utils.hpp>
#include <algorithm>  // std::sort

namespace mesh
{

std::size_t MAX_HASHED_VERTICES = estimate_max_vertices(6);
int INTERNAL_FACE_ID = 0;


std::size_t estimate_max_vertices(const int n_polygon_vertices)
{
  const std::size_t max_vertices =
      static_cast<std::size_t>(
          pow(std::pow(2, 256) - 1, 1.0 / n_polygon_vertices));
  return max_vertices;

}


Point get_element_center(const angem::PointSet<3,double> & vertices,
                         const std::vector<std::size_t>  & ivertices)
{
  return std::move(angem::compute_center_mass
                   (get_vertex_coordinates(vertices, ivertices)));
}


std::vector<Point> get_vertex_coordinates(const angem::PointSet<3,double> & vertices,
                                          const std::vector<std::size_t>  & ivertices)
{
  std::vector<Point> element_vertices(ivertices.size());
  for (std::size_t i=0; i<ivertices.size(); ++i)
    element_vertices[i] = vertices[ivertices[i]];
  return std::move(element_vertices);
}


uint256_t hash_value(const std::vector<std::size_t> & face)
{
  std::vector<std::size_t> face_sorted = face;
  std::sort(face_sorted.begin(), face_sorted.end());

  uint256_t mult = 1;
  uint256_t result = 0;
  for (int i=0; i<face_sorted.size(); ++i)
  {
    // hashing starts with 1 since
    // (1,2,3) and (0,1,2,3) are different elementes
    result += mult * static_cast<uint256_t>(1 + face_sorted[i]);
    mult *= MAX_HASHED_VERTICES;
  }
  return result;
}


int get_face_marker(const uint256_t & hash,
                    const std::unordered_map<uint256_t, int> & map_physical_faces)
{
  const auto it = map_physical_faces.find(hash);
  if (it == map_physical_faces.end())
    return INTERNAL_FACE_ID;
  else return it->second;
}

}
