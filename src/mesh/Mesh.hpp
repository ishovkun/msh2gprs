#pragma once

// internal
#include <mesh_methods.hpp>
#include <ShapeID.hpp>
#include <Face.hpp>
#include <cell_iterator.hpp>
#include <const_cell_iterator.hpp>
#include <face_iterator.hpp>
#include <const_face_iterator.hpp>
// external
#include <angem/Point.hpp>
#include <angem/Polyhedron.hpp>
#include <angem/Polygon.hpp>
// standard
#include <algorithm> // std::sort
#include <unordered_set>


namespace mesh
{

using Point = angem::Point<3,double>;
using Polyhedron = angem::Polyhedron<double>;
using Polygon = angem::Polygon<double>;
using FaceiVertices = std::vector<std::size_t>;


/* This class implements a structure for unstructure grid storage
 * It features constant lookup and insertion times
 * for faces, cells, and their neighbors.
 * Due to this feature though, the order of faces is not guaranteed
 * so that unordered_map[hash]->value should be used to store
 * data related to faces
 */
class Mesh
{
 public:
  Mesh();
  // insert a cell assigned as angem::Polygon
  // Note: this method does not allow duplicates in vertices
  // infers type by number of vertices
  void insert(const Polyhedron & poly,
              const int          marker = -1);
  // insert a cell element assigned as vertex global indices
  void insert_cell(const std::vector<std::size_t> & ivertices,
                   const int                        vtk_id,
                   const int                        marker = -1);
  // insert marker into map_physical_faces
  void insert(const Polygon & poly,
              const int       marker);
  // insert element without searching vertices
  void insert_face(const std::vector<std::size_t> & ivertices,
                   const int                        vtk_id,
                   const int                        marker);

  // ITERATORS

  // Helper function to create cell iterators.
  // Still thinking whether it should be a public method
  cell_iterator create_cell_iterator(const std::size_t icell)
  {return cell_iterator(icell, vertices, cells, map_faces,
                        shape_ids, cell_markers);}
  // create cell iterator for the first cell
  cell_iterator begin_cells(){return create_cell_iterator(0);}
  // end iterator for cells. Must only be used as the range indicator
  cell_iterator end_cells()  {return create_cell_iterator(cells.size());}
  // CONST_ITERATORS
  // Helper function to create cell const_iterators.
  const_cell_iterator create_const_cell_iterator(const std::size_t icell) const
  {return const_cell_iterator(icell, vertices, cells, map_faces,
                              shape_ids, cell_markers);}
  // create cell iterator for the first cell
  const_cell_iterator begin_cells() const {return create_const_cell_iterator(0);}
  // end iterator for cells. Must only be used as the range indicator
  const_cell_iterator end_cells() const {return create_const_cell_iterator(cells.size());}

  // face iterators
  // A helper funciton to create face iterators
 private:
  face_iterator create_face_iterator(const FaceMap::iterator & it)
  {return face_iterator(it, vertices);}
 public:
  // create a face iterator
  face_iterator begin_faces(){return create_face_iterator(map_faces.begin());}
  // create a end iterator for faces
  face_iterator end_faces()  {return create_face_iterator(map_faces.end());}

 private:
  const_face_iterator create_const_face_iterator(FaceMap::const_iterator & it) const
  {return const_face_iterator(it, vertices);}
 public:
  // create a face const_iterator
  const_face_iterator begin_faces() const {return const_face_iterator(map_faces.cbegin(), vertices);}
  // create a end const_iterator for faces
  const_face_iterator end_faces()  const {return const_face_iterator(map_faces.cend(), vertices);}

  // GETTERS
  // get vector of all the grid vertex node coordinates
  std::vector<angem::Point<3,double>> & get_vertices() {return vertices.points;}
  // get const vector of all the grid vertex node coordinates
  const std::vector<angem::Point<3,double>> & get_vertices() const {return vertices.points;}
  // get vertex coordinates
  const angem::Point<3,double> & vertex_coordinates(const std::size_t i) const;
  // get vertex coordinates
  angem::Point<3,double> & vertex_coordinates(const std::size_t i);
  // get vector of neighbor cell indices
  std::vector<std::size_t> get_neighbors( const std::size_t icell ) const;
  // vector of indices of cells neighboring a face
  const std::vector<std::size_t> & get_neighbors( const FaceiVertices & face ) const;
  // get vector of vectors of indices representing faces of a cell
  std::vector<std::vector<std::size_t>> get_faces( const std::size_t ielement ) const;
  // true if vector of cells is empty
  bool empty() const {return cells.empty();}
  // get number of cells
  inline std::size_t n_cells() const {return cells.size();}
  // get number of vertices
  inline std::size_t n_vertices() const {return vertices.size();}
  // get number of faces
  inline std::size_t n_faces() const {return map_faces.size();}
  // get cell center coordinates
  Point get_center(const std::size_t icell) const;
  std::unique_ptr<Polyhedron> get_polyhedron(const std::size_t icell) const;
  // get vector of faces ordered by index (super expernsive -- linear O(n_faces))
  std::vector<face_iterator> get_ordered_faces();

  // MANIPULATION
  // delete an element from the mesh
  void delete_element(const std::size_t ielement);
  // merges jcell into icell if they have a common face
  std::size_t merge_cells(const std::size_t icell,
                          const std::size_t jcell);
  // tell grid which faces to split before calling split_faces method
  void mark_for_split(const face_iterator & face);
  // split faces marked for splitting with mark_for_split
  // returns SurfaceMesh of master DFM faces
  // cleans marked_for_split array upon completion
  SurfaceMesh<double> split_faces();

  // ATTRIBUTES
  angem::PointSet<3,double>             vertices;      // vector of vertex coordinates
  std::vector<std::vector<std::size_t>> cells;         // vertex indices
  std::unordered_map<hash_type, Face>   map_faces;     // map face -> neighbor elements
  std::vector<int>                      shape_ids;     // vector of cell vtk indices
  std::vector<int>                      cell_markers;  // vector of cell markers

 private:
  /* split a vertex
   * retults in adding new vertices (pushed to the back of vertices set)
   * modifies elements of cells array
   */
  void split_vertex(const std::size_t                              ivertex,
                    const std::vector<face_iterator>             & vertex_faces,
                    std::unordered_map<std::size_t, std::size_t> & map_old_new_cells,
                    std::vector<std::vector<std::size_t>>        & new_cells);

  // two elements are in the same group if they are neighbors and
  // the neighboring face is not in vertex_faces array
  std::vector<std::vector<std::size_t>>
  group_cells_based_on_split_faces(const std::unordered_set<std::size_t> & affected_cells,
                                   const std::vector<face_iterator>      & vertex_faces) const;

  // get global indices of polygon face vertices
  std::vector<std::vector<std::size_t>> get_faces(const Polyhedron & poly) const;

  // vector of faces that are markerd for split by the user via mark_for_split
  // Note: the vector is cleared after split_faces is performed
  std::vector<hash_type> marked_for_split;
};


}
