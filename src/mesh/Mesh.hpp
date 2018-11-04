#pragma once

#include <angem/Point.hpp>
#include <angem/Polyhedron.hpp>
#include <uint256/uint256_t.h>

#include <algorithm> // std::sort


namespace mesh
{

using Point = angem::Point<3,double>;
// using Polygon = angem::Polygon<double>;
using Polyhedron = angem::Polyhedron<double>;
using Face = std::vector<std::size_t>;


// maximum 6-vert polygons as faces (hashing)
// this allows hashing about (2^256)^(1/6) ≈ 7x10¹² vertices
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
  std::vector<std::size_t> get_neighbors( const Face & face ) const;
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
  uint256_t hash_value(const Face & face) const;
  const std::size_t max_vertices;
  static constexpr short max_polygon_size_vertices = 6;
};


}
