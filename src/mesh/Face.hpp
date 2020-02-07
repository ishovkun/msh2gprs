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
       const std::vector<std::size_t>        & face_vertices,
       const int                               face_vtk_id,
       const int                               face_marker,
       std::vector<Cell>                     & grid_cells,
       std::vector<Face>                     & grid_faces,
       std::vector<Point>                    & grid_vertices,
       std::vector<std::vector<std::size_t>> & grid_vertex_cells,
       const std::size_t                       parent = constants::default_face_marker);
  // assignment operator
  Face & operator=(const Face & other);
  // comparison operator
  inline bool operator==(const Face & other) const { return index() == other.index(); }
  // comparison operator
  inline bool operator!=(const Face & other) const { return !(*this == other); }
  // ACCESS
  // get face marker
  inline int marker() const { return m_marker; }
  // get face index
  inline std::size_t index() const { return m_index; }
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
  // check if face contains an edge
  bool has_edge(const vertex_pair & edge) const;
  // get the coordinates of face vertices
  std::vector<Point> vertex_coordinates() const;
  // get vector of active neighboring cells.
  std::vector<Cell*> neighbors();
  // get vector of active neighboring cells.
  std::vector<const Cell*> neighbors() const;
  // get vector of neighboring const raw cells.
  std::vector<const Cell*> raw_neighbors() const;
  // get face center of mass
  Point center() const;
  // get face normal
  Point normal() const;
  // return a polygon that forms a face
  Polygon polygon() const;
  // return face polygon area
  double area() const;
  // returns if the face contains a vertex
  // NOTE: complexity is O(n_face_vertices)
  bool has_vertex(const std::size_t vertex_index) const;
  // returns true if has no children; else returns false
  inline bool is_active() const {return m_children.empty();}
  // get reference to parent face. if no parent, simply return itself
  const Face & parent() const { return (*pm_grid_faces)[m_parent]; }
  // get const reference to parent face. if no parent, simply return itself
  Face & parent() { return (*pm_grid_faces)[m_parent]; }
  const Face & ultimate_parent() const;
  Face & ultimate_parent();

  // ---------------------------------- Setters ---------------------------- //

  //    set face marker
  void set_marker(const int marker) { m_marker = marker; }

 protected:
  std::size_t m_index;              // face index
  std::vector<std::size_t> m_vertices;  // face vertex indices
  int m_vtk_id;                         // face vtk id
  int m_marker;                         // face marker
  std::size_t m_parent;                 // parent face
  std::vector<std::size_t> m_children;  // child faces
  // grid stuff
  std::vector<Cell> * pm_grid_cells;  // all grid cells
  std::vector<Face> * pm_grid_faces;                            // grid faces
  std::vector<Point> * pm_grid_vertices;                        // grid vertex coordinates
  std::vector<std::vector<std::size_t>> *pm_grid_vertex_cells;  // vertex -> cell indices
  friend class Mesh;
};

}




