#pragma once

#include <angem/Point.hpp>
#include <angem/Polyhedron.hpp>
#include <uint256/uint256_t.h>
#include <ShapeID.hpp>

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
  // this method does not allow duplicates in vertices
  // infers type by number of vertices
  void insert(const Polyhedron & poly,
              const int          marker = -1);
  bool empty() const {return cells.empty();}

  // GETTERS
  // get vector of neighbor cell indices
  std::vector<std::size_t> get_neighbors( const std::size_t icell ) const;
  // vector of indices of cells neighboring a face
  const std::vector<std::size_t> & get_neighbors( const Face & face ) const;
  // get vector of vectors of indices representing faces of a cell
  std::vector<std::vector<std::size_t>> get_faces( const std::size_t ielement ) const;

  // MANIPULATION
  // delete an element from the mesh
  void delete_element(const std::size_t ielement);
  // merges jcell into icell if they have a common face
  std::size_t merge_cells(const std::size_t icell,
                          const std::size_t jcell);

  // ATTRIBUTES
  angem::PointSet<3,double>             vertices;
  std::vector<std::vector<std::size_t>> cells;  // vertex indices
  // map face -> neighbor elements
  std::unordered_map<uint256_t, std::vector<std::size_t>> map_faces;
  std::vector<int> shape_ids;

 private:
  // get global indices of polygon face vertices
  std::vector<std::vector<std::size_t>> get_faces(const Polyhedron & poly) const;
  // constant complexity (order of n_vertices)
  // linear complexity (size of face)
  uint256_t hash_value(const Face & face) const;
  const std::size_t max_vertices;
  static constexpr int max_polygon_size_vertices = 6;
};


}
