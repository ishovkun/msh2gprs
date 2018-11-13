#pragma once

#include <angem/Point.hpp>
#include <angem/Polyhedron.hpp>
#include <angem/Polygon.hpp>
#include <uint256/uint256_t.h>
#include <ShapeID.hpp>
#include <cell_iterator.hpp>
#include <face_iterator.hpp>
#include <mesh_methods.hpp>

#include <algorithm> // std::sort


/*
 * This class implements a structure for unstructure grid storage
 * It features constant lookup and insertion times
 * for faces, cells, and their neighbors.
 * Due to this feature though, the order of faces is not guaranteed
 * so that unordered_map[hash]->value should be used to store
 * data related to faces
 */

namespace mesh
{

using Point = angem::Point<3,double>;
// using Polygon = angem::Polygon<double>;
using Polyhedron = angem::Polyhedron<double>;
using Polygon = angem::Polygon<double>;
using Face = std::vector<std::size_t>;
using FaceMap = std::unordered_map<uint256_t, std::vector<std::size_t>>;


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
  cell_iterator create_cell_iterator(const std::size_t icell)
  {return cell_iterator(icell, vertices, cells, map_faces,
                        shape_ids, cell_markers);}
  cell_iterator begin_cells(){return create_cell_iterator(0);}
  cell_iterator end_cells()  {return create_cell_iterator(cells.size());}

  // face iterators
  face_iterator create_face_iterator(const FaceMap::iterator & it)
  {return face_iterator(it, vertices, map_physical_faces);}
  face_iterator begin_faces(){return create_face_iterator(map_faces.begin());}
  face_iterator end_faces()  {return create_face_iterator(map_faces.end());}

  // GETTERS
  // get vector of neighbor cell indices
  std::vector<std::size_t> get_neighbors( const std::size_t icell ) const;
  // vector of indices of cells neighboring a face
  const std::vector<std::size_t> & get_neighbors( const Face & face ) const;
  // get vector of vectors of indices representing faces of a cell
  std::vector<std::vector<std::size_t>> get_faces( const std::size_t ielement ) const;
  // true if vector of cells is empty
  bool empty() const {return cells.empty();}
  // get number of cells
  inline std::size_t n_cells() const {return cells.size();}
  // get number of vertices
  inline std::size_t n_vertices() const {return vertices.size();}
  // get cell center coordinates
  Point get_center(const std::size_t icell) const;
  Polyhedron get_polyhedron(const std::size_t icell) const;

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
};


}
