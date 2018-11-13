#pragma once

#include <uint256/uint256_t.h>
#include <angem/Point.hpp>
#include <angem/PointSet.hpp>
#include <Face.hpp>
#include <vector>
#include <unordered_map>

namespace mesh
{
using Point = angem::Point<3,double>;
using FaceMap = std::unordered_map<uint256_t, Face>;

class face_iterator
{
 public:
  face_iterator(const FaceMap::iterator            & it,
                angem::PointSet<3,double>          & vertices);
  // comparison
  bool operator==(const face_iterator & other) const;
  bool operator!=(const face_iterator & other) const;
  // Getters
  // get center of mass of a face
  Point center() const;
  // get face normal
  Point normal() const;
  // get face marker, (-1) if not defined
  int marker() const;
  // get hash value
  inline uint256_t hash() const {return face_it->first;}
  // get vector of neighbor indices
  inline
  const std::vector<std::size_t> & neighbors() const {return face_it->second.neighbors;}
  std::vector<Point> vertices() const;
  std::vector<std::size_t> vertex_indices() const;
  // incrementing
  face_iterator & operator++();
  face_iterator & operator--();

 private:
  FaceMap::iterator face_it;
  angem::PointSet<3,double>          & mesh_vertices;
};

}
