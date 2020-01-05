#pragma once

#include "constants.hpp"
#include "angem/Point.hpp"
#include "angem/Polygon.hpp"
#include <vector>

namespace mesh
{
using Point = angem::Point<3,double>;
using Polygon = angem::Polygon<double>;
using vertex_pair = std::pair<size_t, size_t>;
class Cell;

class Face
{
 public:
  Face(const std::size_t                       face_index,
       const std::size_t                       master_face_index,
       const std::vector<std::size_t>        & face_vertices,
       const int                               face_vtk_id,
       const int                               face_marker,
       std::vector<Cell>                     & grid_cells,
       std::vector<Point>                    & grid_vertices,
       std::vector<std::vector<std::size_t>> & grid_vertex_cells,
       const std::size_t                       parent = constants::default_face_marker);
  // assignment operator
  Face & operator=(const Face & other);
  // comparison operator
  inline bool operator==(const Face & other) const { return index() == other.index(); }
  // ACCESS
  // get face marker
  inline int marker() const { return m_marker; }
  // get face index
  inline std::size_t index() const { return m_index; }
  // get index of master face (before splitting)
  inline std::size_t master_index() const { return m_master_face_index; }
  // get index of parent face. if no parent, simply return the index
  inline std::size_t parent() const { return m_parent; }
  // get vtk index of a face polygon
  inline int vtk_id() const { return m_vtk_id; }
  // return indices of child faces
  inline const std::vector<std::size_t> & children() const { return m_children; }
  // get the indices of face vertices
  inline std::vector<std::size_t> & vertices() { return m_vertices; }
  // get the indices of face vertices
  inline const std::vector<std::size_t> & vertices() const { return m_vertices; }
  // get vector of index pairs that represent edges
  std::vector<vertex_pair> edges() const;
  // get the coordinates of face vertices
  std::vector<Point> vertex_coordinates() const;
  // get vector of neighboring cells
  std::vector<Cell*> neighbors();
  // get vector of neighboring const cells
  std::vector<const Cell*> neighbors() const;
  // get face center of mass
  Point center() const;
  // get face normal
  Point normal() const;
  Polygon polygon() const;
  double area() const;
  bool has_vertex(const std::size_t vertex_index) const;
  // returns true if has no children; else returns false
  inline bool is_active() const {return m_children.empty();}

  // ---------------------------------- Setters ---------------------------- //

  //  set face marker
  void set_marker(const int marker) { m_marker = marker; }

 protected:
  std::size_t m_index;              // face index
  std::size_t m_master_face_index;  // index before split
  std::vector<std::size_t> m_vertices;  // face vertex indices
  int m_vtk_id;                         // face vtk id
  int m_marker;                         // face marker
  std::size_t m_parent;                 // parent face
  std::vector<std::size_t> m_children;  // child faces
  // grid stuff
  std::vector<Cell> * pm_grid_cells;  // all grid cells
  std::vector<Point> * pm_grid_vertices;                        // grid vertex coordinates
  std::vector<std::vector<std::size_t>> *pm_grid_vertex_cells;  // vertex -> cell indices
  friend class Mesh;
};

}




