#pragma once

#include <uint256/uint256_t.h>
#include <angem/PointSet.hpp>
#include <Polyhedron.hpp>
#include <Face.hpp>
#include <face_iterator.hpp>
#include <vector>
#include <unordered_map>

namespace mesh
{
using Point = angem::Point<3,double>;

class cell_iterator
{
 public:
  // constructor
  cell_iterator(const std::size_t                       icell,
                angem::PointSet<3,double>             & vertices,
                std::vector<std::vector<std::size_t>> & cells,
                std::unordered_map<uint256_t, Face>   & map_faces,
                std::vector<int>                      & shape_ids,
                std::vector<int>                      & cell_markers);
  // assignment operator
  cell_iterator & operator=(const cell_iterator & other);
  // comparison
  bool operator==(const cell_iterator & other) const;
  bool operator!=(const cell_iterator & other) const;
  // GETTERS
  inline std::size_t index() const {return icell;}
  inline int shape_id() const {return shape_ids[icell];}
  inline int marker() const {return cell_markers[icell];}
  std::vector<std::size_t> neighbor_indices() const;
  std::vector<cell_iterator> neighbors() const;
  Point center() const;
  angem::Polyhedron<double> polyhedron() const;
  std::vector<face_iterator> faces() const;


  // incrementing
  cell_iterator & operator++();
  cell_iterator & operator--();

 private:
  // std::vector<std::vector<std::size_t>> & cells;  // vertex indices
  std::size_t icell;
  angem::PointSet<3,double>             & vertices;
  std::vector<std::vector<std::size_t>> & cells;
  std::unordered_map<uint256_t, Face>   & map_faces;
  std::vector<int> & shape_ids;
  std::vector<int> & cell_markers;
};

}
