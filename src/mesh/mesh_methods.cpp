#include <mesh_methods.hpp>
#include <angem/utils.hpp>
#include <algorithm>  // std::sort

namespace mesh
{

const int MAX_POLYGON_VETRICES = 6;
const std::size_t MAX_HASHED_VERTICES = estimate_max_vertices(MAX_POLYGON_VETRICES);
const int INTERNAL_FACE_ID = 0;


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
  for (int i=face_sorted.size()-1; i>=0; --i)
  {
    // hashing starts with 1 since
    // (1,2,3) and (0,1,2,3) are different elementes
    result += mult * static_cast<uint256_t>(1 + face_sorted[i]);
    mult *= MAX_HASHED_VERTICES;
  }
  return result;
}


std::vector<std::size_t> invert_hash(const uint256_t & hash)
{
  std::vector<std::size_t> vertices;
  uint256_t current_value = hash;
  for (int i=0; i<MAX_POLYGON_VETRICES; ++i)
  {
    // - 1 since hashing starts with 1
    const std::size_t ivertex = (current_value % MAX_HASHED_VERTICES) - 1;
    current_value = current_value / MAX_HASHED_VERTICES;
    vertices.push_back(ivertex);
    if (current_value == 0)
      break;
  }
  std::vector<std::size_t> vertices_sorted(vertices.size());
  int ind = vertices.size() - 1;
  for (const auto & v : vertices)
  {
    vertices_sorted[ind] = v;
    ind--;
  }
  return std::move(vertices_sorted);
}


int get_face_marker(const uint256_t & hash,
                    const std::unordered_map<uint256_t, int> & map_physical_faces)
{
  const auto it = map_physical_faces.find(hash);
  if (it == map_physical_faces.end())
    return INTERNAL_FACE_ID;
  else return it->second;
}


std::vector<std::vector<std::size_t>>
get_face_indices(const angem::Polyhedron<double> & poly,
                 const angem::PointSet<3,double> & vertices)
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


// angem::Polyhedron<double> get_polyhedron(const std::size_t icell,
//                                          const angem::PointSet<3,double> & vertices,
//                                          const std::vector<int>          & shape_ids)
// {
//   return std::move(angem::PolyhedronFactory::create<double>(vertices.points,
//                                                             cells[icell],
//                                                             shape_ids[icell]));
// }

}
