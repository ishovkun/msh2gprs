#pragma once

#include "angem/Point.hpp"
#include "angem/Polyhedron.hpp"
#include "Face.hpp"
#include <vector>
#include <functional>  // std::reference_wrapper

namespace mesh
{
using Point = angem::Point<3,double>;
using Polyhedron = angem::Polyhedron<double>;

static const int DEFAULT_CELL_MARKER = -1;

class Cell
{
 public:
  // constructor
  Cell(const std::size_t          cell_index,
       const std::vector<std::size_t> & vertices,
       const std::vector<std::size_t> & faces,
       std::vector<Point>        & grid_vertices,
       std::vector<Cell>        & grid_cells,
       std::vector<Face>        & grid_faces,
       const int                  vtk_id,
       const int                  marker);
  // access operators
  // get cell index
  inline std::size_t index() const { return m_index; }
  // get cell marker
  inline int marker() const { return m_marker; }
  // get vtk id
  inline int vtk_id() const { return m_vtk_id; }
  // get vector of vertex indices
  inline std::vector<std::size_t> & vertices() { return m_vertices; }
  // get const vector of vertex indices
  inline const std::vector<std::size_t> & vertices() const { return m_vertices; }
  // get vector of neighbors
  std::vector<Cell*> neighbors();
  // get vector of neighbors
  std::vector<const Cell*> neighbors() const;
  // get vector of cell faces
  std::vector<Face*> faces();
  // get vector of cell faces
  std::vector<const Face*> faces() const;
  // convenience methods
  // get the number of cell vertices
  inline std::size_t n_vertices() const { return m_vertices.size(); }
  // get center of mass
  Point center() const;
  // get cell volume
  double volume() const;
  // get a polyhedron that represents a cell
  std::unique_ptr<Polyhedron> polyhedron() const;
  // true if cell hace a vertex, false otherwise
  bool has_vertex(const std::size_t vertex_index) const;

 protected:
  // this cell stuff
  std::size_t m_index;
  int m_marker, m_vtk_id;
  std::vector<std::size_t> m_vertices;
  std::vector<std::size_t> m_faces;
  // grid stuff
  std::vector<Point> & m_grid_vertices;
  std::vector<Cell> & m_grid_cells;
  std::vector<Face> & m_grid_faces;
  friend class Mesh;
};

}  // end namespace