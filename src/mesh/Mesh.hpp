#pragma once

#include <angem/Point.hpp>
#include <angem/Polygon.hpp>

#include <algorithm> // std::sort

using Point = angem::Point<3,double>;
using Polygon = angem::Polygon<double>;
using Polyhedron = angem::Polyhedron<double>;
using Face = std::vector<std::size_t>;


namespace mesh
{

using namespace angem;

// N is the maximum number of face polygon vertices
// it determines the hashing algorithm
// the less N is, more faces can be hashed
// 4 is the recommended number
template <short N>
class Mesh
{
 public:
  Mesh();
  void insert(const Polyhedron & poly);
  bool empty() const {return cells.empty();}

  // GETTERS
  // get vector of neighbor cell indices
  std::vector<std::size_t> get_neighbors( const std::size_t icell ) const;
  // vector of indices of cells neighboring a face
  std::vector<std::size_t> get_neighbors( const Polygon & edge ) const;
  // get vector of vectors of indices representing faces of a cell
  std::vector<std::vector<std::size_t>> get_faces( const std::size_t ielement ) const;

  // MANIPULATION
  // delete an element from the mesh
  void delete_element(const std::size_t ielement);
  // merges jcell into icell if they have a common face
  std::size_t merge_cells(const std::size_t icell,
                          const std::size_t jcell);

  // ATTRIBUTES
  angem::PointSet<3,double>                    vertices;
  std::vector<std::vector<std::size_t>> cells;  // indices
  // hash of two vert indices -> vector polygons
  // essentially edge -> neighbor elements
  std::unordered_map<std::size_t, std::vector<std::size_t>> map_faces;

 private:
  // constant complexity (order of n_vertices)
  // linear complexity (size of face)
  __int128 hash_value(const Face & face) const;
  const std::size_t max_vertices =
      static_cast<std::size_t>(
          pow(std::numeric_limits<__int128>::max() - 1, 1.0 / N));
};


template <short N>
Mesh<N>::Mesh()
{
  static_assert(N > 0);
}


template <short N>
__int128 Mesh<N>::hash_value(const Face & face) const
{
  assert(face.size() < N);

  Face face_sorted = face;
  std::sort(face_sorted.begin(), face_sorted.end());

  __int128 mult = 1;
  __int128 result = 0;
  for (short i=0; i<N; ++i)
  {
    // hashing starts with 1 since
    // (1,2,3) and (0,1,2,3) are different elementes
    result += 1 + mult * face_sorted[i];
    mult *= max_vertices;
  }
  return result;
}

template <short N>
void Mesh<N>::insert(const Polyhedron & poly)
{
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
