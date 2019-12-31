#pragma once

#include "Face.hpp"
#include "Cell.hpp"
#include "active_cell_iterator.hpp"
#include "active_cell_const_iterator.hpp"
#include "SurfaceMesh.hpp"
#include <angem/Point.hpp>
#include <angem/Polyhedron.hpp>
#include <angem/Polygon.hpp>
// #include <cell_const_iterator.hpp>
// #include <face_iterator.hpp>
// #include <face_const_iterator.hpp>
#include <algorithm> // std::sort
#include <unordered_set>


namespace mesh
{

static const int DEFAULT_FACE_MARKER = -1;
using Point = angem::Point<3,double>;
using Polyhedron = angem::Polyhedron<double>;
using Polygon = angem::Polygon<double>;
using FaceiVertices = std::vector<std::size_t>;

/* This class implements a structure for unstructure grid storage
 * It features constant lookup and insertion times
 * for faces, cells, and their neighbors.
 * Due to this feature though, the order of faces is not guaranteed
 * so that unordered_map[hash]->value should be used to store
 * data related to faces */
class Mesh
{
 public:
  Mesh();
  /* Insert a standard vtk polyhedron cell into grid */
  std::size_t insert_cell(const std::vector<std::size_t> & ivertices,
                          const int                        vtk_id,
                          const int                        marker = DEFAULT_CELL_MARKER);
  /* insert a face into the grid */
  std::size_t insert_face(const std::vector<std::size_t> & ivertices,
                          const int                        vtk_id,
                          const int                        marker = DEFAULT_CELL_MARKER);

  // ITERATORS
  //  create cell iterator for the first active cell
  active_cell_iterator begin_active_cells();
  //  create cell iterator for the first active cell
  active_cell_const_iterator begin_active_cells() const;
  // Returns an iterator referring to the past-the-end active cell
  inline active_cell_iterator end_active_cells() {return active_cell_iterator(nullptr);}
  // Returns an iterator referring to the past-the-end active cell
  inline active_cell_const_iterator end_active_cells() const {return active_cell_const_iterator(nullptr);}
  /* RAW cell iterators, use with caution */
  //  create cell iterator for the first cell
  //  NOTE: RAW cell iterator, use with caution
  inline std::vector<Cell>::iterator begin_cells() {return m_cells.begin();}
  // create cell const_iterator for the first cell
  //  NOTE: RAW cell iterator, use with caution
  inline std::vector<Cell>::const_iterator begin_cells() const {return m_cells.begin();}
  // Returns an iterator referring to the past-the-end element in the cell container
  //  NOTE: RAW cell iterator, use with caution
  inline std::vector<Cell>::iterator end_cells()   {return m_cells.end();}
  // end iterator for cells. Must only be used as the range indicator
  //  NOTE: RAW cell iterator, use with caution
  inline std::vector<Cell>::const_iterator end_cells() const {return m_cells.end();}
  // create a face iterator
  inline std::vector<Face>::iterator begin_faces(){return m_faces.begin();}
  // create a face const_iterator
  inline std::vector<Face>::const_iterator begin_faces() const {return m_faces.begin();}
  // create an end iterator for faces
  inline std::vector<Face>::iterator end_faces(){return m_faces.end();}
  // create an end const_iterator for faces
  inline std::vector<Face>::const_iterator end_faces() const {return m_faces.end();}
  // create a vertex iterator pointing to the first vertex
  inline std::vector<Point>::iterator begin_vertices() { return m_vertices.begin(); }
  // create a vertex const_iterator pointing to the first vertex
  inline std::vector<Point>::const_iterator begin_vertices() const { return m_vertices.begin(); }
  // create a vertex iterator pointing to the past last vertex
  inline std::vector<Point>::iterator end_vertices() { return m_vertices.end(); }
  // create a vertex const_iterator pointing to the past last vertex
  inline std::vector<Point>::const_iterator end_vertices() const { return m_vertices.end(); }

  // ACCESS OPERATORS

  // get reference to a cell object by index
  inline Cell & cell(const size_t cell_index) { return m_cells[cell_index]; }
  // get const reference to a cell object by index
  inline const Cell & cell(const size_t cell_index) const { return m_cells[cell_index]; }
  // get reference to a face object by index
  inline Face & face(const size_t face_index) { return m_faces[face_index]; }
  // get const reference to a face object by index
  inline const Face & face(const size_t face_index) const { return m_faces[face_index]; }
  // get reference to a vertex object by index
  inline Point& vertex(const size_t vertex_index) {return m_vertices[vertex_index];}
  // get const reference to a vertex object by index
  inline const Point& vertex(const size_t vertex_index) const {return m_vertices[vertex_index];}
  // GETTERS
  // get vector of all the grid vertex node coordinates
  inline std::vector<angem::Point<3,double>> & vertices() {return m_vertices;}
  // get const vector of all the grid vertex node coordinates
  inline const std::vector<angem::Point<3,double>> & vertices() const {return m_vertices;}
  // get vector of cells
  inline std::vector<Cell> & cells() {return m_cells;}
  // get vector of cells
  inline const std::vector<Cell> & cells() const {return m_cells;}
  // true if vector of cells is empty
  bool empty() const {return m_cells.empty();}
  // get number of cells
  inline std::size_t n_cells() const {return m_cells.size();}
  // get number of vertices
  inline std::size_t n_vertices() const {return m_vertices.size();}
  // get number of faces
  inline std::size_t n_faces() const {return m_faces.size();}

  // MANIPULATION
  // delete a cell mesh
   void delete_cell(const std::size_t ielement);
  //  // merges jcell into icell if they have a common face
  std::size_t merge_cells(const std::size_t icell, const std::size_t jcell);
  // tell grid which faces to split before calling split_faces method
  void mark_for_split(const std::size_t face_index);
  // split faces marked for splitting with mark_for_split
  // returns SurfaceMesh of master DFM faces
  // cleans marked_for_split array upon completion
  SurfaceMesh<double> split_faces();

  /* Split a cell by cutting it with a plane. New cell indices are appended
   * at the back, so it is safe to split multiple cells in a row. */
  void split_cell(Cell & cell, const angem::Plane<double> & plane);
  /* split a vertex
   * retults in adding new vertices (pushed to the back of vertices set) */
  void split_vertex(const std::size_t               vertex_index,
                    const std::vector<std::size_t> &splitted_face_indices);

 private:
  void insert_vertex(const std::size_t vertex,
                     const std::size_t face,
                     const std::size_t cell);

  // two elements are in the same group if they are neighbors and
  // the neighboring face is not in vertex_faces array
  std::vector<std::vector<std::size_t>>
  group_cells_based_on_split_faces(const std::vector<std::size_t> & affected_cells,
                                   const std::vector<std::size_t> & splitted_face_indices) const;


  /* private insert cell class that does all the cell insertion work */
  std::size_t insert_cell_(const std::vector<std::size_t> & ivertices,
                           const std::vector<std::vector<std::size_t>> & cell_faces,
                           const int                        vtk_id,
                           const int                        marker);
  /* Insert an arbitrary polyhedron cell into grid .
   * A wrapper on the above function to minimize bookkeeping. */
  std::size_t insert_cell_(const std::vector<std::vector<std::size_t>> & cell_faces,
                           const int                        marker = DEFAULT_CELL_MARKER);


 private:
  // ATTRIBUTES
  std::vector<angem::Point<3,double>>   m_vertices;      // vector of vertex coordinates
  std::vector<Cell>                     m_cells;
  std::vector<Face>                     m_faces;
  std::vector<std::vector<std::size_t>> m_vertex_cells;  // vertex neighboring cells
  std::vector<std::vector<std::size_t>> m_vertex_faces;  // vertex neighboring faces
  // Used as a tmp container when splitting faces for dfm
  std::vector<std::size_t> m_faces_marked_for_split;
};


}
