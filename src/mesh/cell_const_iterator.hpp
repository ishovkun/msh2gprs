#pragma once

#include "angem/Polyhedron.hpp"
#include "Cell.hpp"
#include "Face.hpp"
#include "face_const_iterator.hpp"
#include <vector>
#include <unordered_map>
#include <mesh_methods.hpp>

#include <memory>  // unique_ptr

namespace mesh
{
using Point = angem::Point<3,double>;
using Polyhedron = angem::Polyhedron<double>;

class cell_const_iterator
{
 public:
  // constructor
  cell_const_iterator(const std::size_t                           cell_index,
                      const std::vector<Cell>                   & grid_cells,
                      const std::vector<Face>                   & grid_faces,
                      const std::vector<angem::Point<3,double>> & grid_vertices);
  // access operator
  inline const Cell * operator*() const {
    if (cell_index >= grid_cells.size())
      throw std::out_of_range("cell iterator past end");
    return &(grid_cells[cell_index]);
  }

  // assignment operator
  cell_const_iterator & operator=(const cell_const_iterator & other);
  // comparison
  bool operator==(const cell_const_iterator & other) const;
  bool operator!=(const cell_const_iterator & other) const;

  // GETTERS
  inline const std::vector<std::size_t> & vertex_indices() const {return operator*()->vertices;}
  inline size_t n_vertices() const {return vertex_indices().size();}
  cell_const_iterator neighbor_by_face(const face_const_iterator & face) const;
  inline std::size_t index() const {return cell_index;}
  inline int vtk_id() const {return operator*()->vtk_id;}
  inline int marker() const {return operator*()->marker;}
  std::vector<std::size_t> neighbor_indices() const;
  std::vector<cell_const_iterator> neighbors() const;
  Point center() const;
  double volume() const;
  std::unique_ptr<Polyhedron> polyhedron() const;
  const std::vector<Face*> faces() const;

  // other
  bool has_vertex(const std::size_t ivertex) const;

  // incrementing
  cell_const_iterator & operator++();
  cell_const_iterator & operator--();

 private:
  std::size_t                                   cell_index;
  const std::vector<angem::Point<3,double>>   & grid_vertices;
  const std::vector<Cell>                     & grid_cells;
  const std::vector<Face>                     & grid_faces;
};

} // end namespace
