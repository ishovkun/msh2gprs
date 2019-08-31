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

class face_const_iterator : public std::iterator<std::random_access_iterator_tag, Face>
{
 public:
  // Default constructor
  // Default constructor
  face_const_iterator(const std::size_t                           face_index,
                      const std::vector<Face>                   & grid_faces,
                      const std::vector<angem::Point<3,double>> & grid_vertices);
  // Copy constructor
  face_const_iterator(const face_const_iterator & other);

  // access operator
  const Face * operator*() const { return &(grid_faces[face_index]); }
  // comparison
  // returns true if both iterators point to the same cell
  bool operator==(const face_const_iterator & other) const;
  // returns true if the iterators point to different cells
  bool operator!=(const face_const_iterator & other) const;

  // SETTERS
  // assignment operator
  face_const_iterator & operator=(const face_const_iterator & other);

  // GETTERS
  // get center of mass of a face
  Point center() const;
  // get face normal vector
  Point normal() const;
  // get face marker, (-1) if not defined
  int marker() const;
  // get face index
  std::size_t index() const {return operator*()->index;}
  // get an index of the parent (master) face that existed before the split
  // same as index() if the face has not been split
  std::size_t master_index() const {return operator*()->master_face_index;}
  // get vtk id of the face
  int vtk_id() const {return operator*()->vtk_id;}
  // get vector of neighbor cell indices
  inline const std::vector<std::size_t> & neighbors() const {return operator*()->neighbor_cells;}
  // get vector of face vertex coordinates
  std::vector<Point> vertex_coordinates() const;
  // get vector of face vertex indices
  std::vector<std::size_t> vertex_indices() const;
  // get angem::Polygon from vertices
  angem::Polygon<double> polygon() const;
  // incrementing
  // increment operator
  face_const_iterator & operator++();
  // decrement operator
  face_const_iterator & operator--();

 private:
  std::size_t                             face_index;
  const std::vector<Face>                     & grid_faces;
  const std::vector<angem::Point<3,double>>   & grid_vertices;
};

}
