#pragma once

#include <angem/PointSet.hpp>
#include <Polyhedron.hpp>
#include <Face.hpp>
#include <const_face_iterator.hpp>
#include <vector>
#include <unordered_map>
#include <mesh_methods.hpp>

#include <memory>  // unique_ptr

namespace mesh
{
using Point = angem::Point<3,double>;

class const_cell_iterator
{
 public:
  // constructor
  const_cell_iterator(const std::size_t                       icell,
                      const angem::PointSet<3,double>             & vertices,
                      const std::vector<std::vector<std::size_t>> & cells,
                      const std::unordered_map<hash_type, Face>   & map_faces,
                      const std::vector<int>                      & shape_ids,
                      const std::vector<int>                      & cell_markers);
  // assignment operator
  const_cell_iterator & operator=(const const_cell_iterator & other);
  // comparison
  bool operator==(const const_cell_iterator & other) const;
  bool operator!=(const const_cell_iterator & other) const;
  // GETTERS
  inline const std::vector<std::size_t> & vertices() {return cells[icell];}
  const_cell_iterator neighbor_by_face(const const_face_iterator & face) const;
  inline std::size_t index() const {return icell;}
  inline int shape_id() const {return shape_ids[icell];}
  inline int marker() const {return cell_markers[icell];}
  std::vector<std::size_t> neighbor_indices() const;
  std::vector<const_cell_iterator> neighbors() const;
  Point center() const;
  double volume() const;
  std::unique_ptr<Polyhedron<double>> polyhedron() const;
  std::vector<const_face_iterator> faces() const;

  // other
  bool has_vertex(const std::size_t ivertex) const;

  // incrementing
  const_cell_iterator & operator++();
  const_cell_iterator & operator--();

 private:
  // std::vector<std::vector<std::size_t>> & cells;  // vertex indices
  std::size_t icell;
  const angem::PointSet<3,double>             & mesh_vertices;
  const std::vector<std::vector<std::size_t>> & cells;
  const std::unordered_map<hash_type, Face>   & map_faces;
  const std::vector<int> & shape_ids;
  const std::vector<int> & cell_markers;
};

} // end namespace
