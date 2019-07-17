#include <mesh_methods.hpp>
#include <angem/utils.hpp>
#include <algorithm>  // std::sort
#include <cmath>      // std::pow

namespace mesh
{

// maximum 6-vert polygons as faces (hashing)
// this allows hashing about (2^256)^(1/6) ≈ 7x10¹² vertices
const int MAX_POLYGON_VETRICES = 6;
const std::size_t MAX_HASHED_VERTICES = estimate_max_vertices(MAX_POLYGON_VETRICES);
const int INTERNAL_FACE_ID = 0;


std::size_t estimate_max_vertices(const int n_polygon_vertices)
{
  const std::size_t max_vertices =
      static_cast<std::size_t>(
          std::pow(std::pow(2, 256) - 1, 1.0 / n_polygon_vertices));
  return max_vertices;

}


Point get_element_center(const angem::PointSet<3,double> & vertices,
                         const std::vector<std::size_t>  & ivertices)
{
  return angem::compute_center_mass(get_vertex_coordinates(vertices, ivertices));
}


std::vector<Point> get_vertex_coordinates(const angem::PointSet<3,double> & vertices,
                                          const std::vector<std::size_t>  & ivertices)
{
  std::vector<Point> element_vertices(ivertices.size());
  for (std::size_t i=0; i<ivertices.size(); ++i)
    element_vertices[i] = vertices[ivertices[i]];
  return element_vertices;
}


std::vector<Point> get_vertex_coordinates(const angem::PointSet<3,double> * vertices,
                                          const std::vector<std::size_t>  & ivertices)
{
  std::vector<Point> element_vertices(ivertices.size());
  for (std::size_t i=0; i<ivertices.size(); ++i)
    element_vertices[i] = (*vertices)[ivertices[i]];
  return element_vertices;
}


hash_type hash_value(const std::vector<std::size_t> & face)
{
  std::vector<std::size_t> face_sorted = face;
  std::sort(face_sorted.begin(), face_sorted.end());

  hash_type mult(1);
  hash_type result(0);
  for (int i=face_sorted.size()-1; i>=0; --i)
  {
    // hashing starts with 1 since
    // (1,2,3) and (0,1,2,3) are different elementes
    result += mult * static_cast<hash_type>(1 + face_sorted[i]);
    mult *= MAX_HASHED_VERTICES;
  }
  return result;
}


std::vector<std::size_t> invert_hash(const hash_type & hash)
{
  std::vector<std::size_t> vertices;
  hash_type current_value = hash;
  hash_type mvert(MAX_HASHED_VERTICES);
  hash_type zero(0);
  for (int i=0; i<MAX_POLYGON_VETRICES; ++i)
  {
    // - 1 since hashing starts with 1
    // const std::size_t ivertex = (current_value % mvert) - 1;
    // [0] i guess that's how to retreive size_t from this class
    const std::size_t ivertex = static_cast<std::size_t>(current_value % mvert) - 1;
    // current_value = current_value / MAX_HASHED_VERTICES;
    current_value = current_value / mvert;
    vertices.push_back(ivertex);
    if (current_value == zero)
      break;
  }

  std::vector<std::size_t> vertices_sorted(vertices.size());
  int ind = vertices.size() - 1;
  for (const auto & v : vertices)
  {
    vertices_sorted[ind] = v;
    ind--;
  }
  return vertices_sorted;
}


int get_face_marker(const hash_type & hash,
                    const std::unordered_map<hash_type, int> & map_physical_faces)
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
  return faces;

}


std::unordered_set<std::size_t> find_internal_vertices(const SurfaceMesh<double> & msh)
{
  std::unordered_set<std::size_t> internal_points;

  for (const_edge_iterator edge=msh.begin_edges(); edge!=msh.end_edges(); ++edge)
  {
    if (edge.neighbors().size() > 1)
    {
      const std::pair<std::size_t, std::size_t> ivertices = edge.vertex_indices();
      std::cout << "edge.neighbors().size() = " << edge.neighbors().size() << std::endl;
      // TODO: record number of neighbors
      // if > 2 -- fracture intersection. need more vertex splits
      internal_points.insert(ivertices.first);
      internal_points.insert(ivertices.second);
    }
  }
  return internal_points;
}

}
