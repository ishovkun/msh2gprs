#pragma once

#include "angem/Polyhedron.hpp"
#include <Face.hpp>
#include <face_iterator.hpp>
#include <vector>
// #include <unordered_map>
// #include <mesh_methods.hpp>
#include "Cell.hpp"

#include <memory>  // unique_ptr

namespace mesh
{
using Point = angem::Point<3,double>;
using Polyhedron = angem::Polyhedron<double>;

class cell_iterator : public std::iterator<std::random_access_iterator_tag, Cell>
{
 public:
  // constructor
  cell_iterator(const std::size_t                     cell_index,
                std::vector<Cell> &                   grid_cells,
                std::vector<Face> &                   grid_faces,
                std::vector<angem::Point<3,double>> & grid_vertices);
  // assignment operator
  cell_iterator & operator=(const cell_iterator & other);
  // comparison
  bool operator==(const cell_iterator & other) const;
  bool operator!=(const cell_iterator & other) const;
  // GETTERS
  // access operator
  inline Cell * operator*() { return &(grid_cells[cell_index]); }
  inline std::size_t n_vertices() { return operator*()->vertices.size(); }
  inline std::vector<std::size_t> & vertex_indices() {return operator*()->vertices;}
  cell_iterator neighbor_by_face(const face_iterator & face) const;
  inline std::size_t index() const {return cell_index;}
  inline int vtk_id() const {return operator*()->vtk_id;}
  inline int marker() const {return operator*()->marker;}
  std::vector<Face*> faces();
  // std::vector<std::size_t> neighbor_indices() const;
  // std::vector<cell_iterator> neighbors() const;
  Point center() const;
  double volume() const;
  std::unique_ptr<Polyhedron> polyhedron() const;

  // other
  bool has_vertex(const std::size_t ivertex) const;

  // incrementing
  cell_iterator & operator++();
  cell_iterator & operator--();

 private:
  inline const Cell * operator*() const { return &(grid_cells[cell_index]); }

  // ATTRIBUTES
 private:
  std::size_t cell_index;
  std::vector<angem::Point<3,double>>   & grid_vertices;
  std::vector<Cell>                     & grid_cells;
  std::vector<Face>                     & grid_faces;
};

}
