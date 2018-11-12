#pragma once

#include <angem/Point.hpp>
#include <angem/Polyhedron.hpp>
#include <angem/Polygon.hpp>
#include <uint256/uint256_t.h>
#include <ShapeID.hpp>
#include <CellIterator.hpp>
#include <FaceIterator.hpp>

#include <algorithm> // std::sort


namespace mesh
{

using Point = angem::Point<3,double>;
// using Polygon = angem::Polygon<double>;
using Polyhedron = angem::Polyhedron<double>;
using Polygon = angem::Polygon<double>;
using Face = std::vector<std::size_t>;

/* TODO: iterators of cells, faces, physical faces
 * functions for element volume and face area and normals
 */

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
  // insert marker into map_physical_faces
  void insert_physical_face(const Polygon & poly,
                            const int       marker);

  // iterators
  // cell iterators
  CellIterator create_cell_iterator(const std::size_t icell)
  {return CellIterator(icell, map_faces, shape_ids, cell_markers);}
  CellIterator begin_cells(){return create_cell_iterator(0);}
  CellIterator end_cells()  {return create_cell_iterator(cells.size());}

  // face iterators
  FaceIterator create_face_iterator(const std::size_t iface)
  {return FaceIterator(iface, map_faces, map_physical_faces);}
  FaceIterator begin_faces(){return create_face_iterator(0);}
  FaceIterator end_faces()  {return create_face_iterator(map_faces.size());}

  // GETTERS
  // get vector of neighbor cell indices
  std::vector<std::size_t> get_neighbors( const std::size_t icell ) const;
  // vector of indices of cells neighboring a face
  const std::vector<std::size_t> & get_neighbors( const Face & face ) const;
  // get vector of vectors of indices representing faces of a cell
  std::vector<std::vector<std::size_t>> get_faces( const std::size_t ielement ) const;
  // true if vector of cells is empty
  bool empty() const {return cells.empty();}
  std::size_t n_cells() const {return cells.size();}
  std::size_t n_vertices() const {return vertices.size();}

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
  // map face -> marker
  std::unordered_map<uint256_t, int> map_physical_faces;
  std::vector<int> shape_ids;
  std::vector<int> cell_markers;

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
