#pragma once

#include <uint256/uint256_t.h>
#include <angem/Point.hpp>
#include <angem/PointSet.hpp>
#include <Face.hpp>
#include <mesh_methods.hpp>
#include <unordered_map>
#include <vector>

namespace mesh
{
using Point = angem::Point<3,double>;
using FaceMap = std::unordered_map<hash_type, Face>;

class face_iterator
{
 public:
  // Default constructor
  face_iterator(const FaceMap::iterator            & it,
                angem::PointSet<3,double>          & vertices);
  // Copy constructor
  face_iterator(const face_iterator & other);

  // comparison
  // returns true if both iterators point to the same cell
  bool operator==(const face_iterator & other) const;
  // returns true if the iterators point to different cells
  bool operator!=(const face_iterator & other) const;

  // SETTERS
  // assignment operator
  face_iterator & operator=(const face_iterator & other);

  // GETTERS
  // get center of mass of a face
  Point center() const;
  // get face normal vector
  Point normal() const;
  // get face marker, (-1) if not defined
  int marker() const;
  // get hash value of the face
  inline hash_type hash() const {return face_it->first;}
  // get face index
  std::size_t index() const {return face_it->second.index;}
  // get an index of the parent (master) face that existed before the split
  // same as index() if the face has not been split
  std::size_t master_index() const {return face_it->second.old_index;}
  // get vtk id of the face
  int vtk_id() const {return face_it->second.vtk_id;}
  // get vector of neighbor cell indices
  inline
  const std::vector<std::size_t> & neighbors() const {return face_it->second.neighbors;}
  // get vector of face vertex coordinates
  std::vector<Point> vertices() const;
  // get vector of face vertex indices
  std::vector<std::size_t> vertex_indices() const;
  // get angem::Polygon from vertices
  angem::Polygon<double> polygon() const;
  // incrementing
  // increment operator
  face_iterator & operator++();
  // decrement operator
  face_iterator & operator--();

 private:
  FaceMap::iterator           face_it;        // iterator in the face container
  angem::PointSet<3,double> & mesh_vertices;  // reference to the vertices container
};

}
