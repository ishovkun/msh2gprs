#pragma once

#include "Face.hpp"

namespace mesh {

class active_face_const_iterator
{
 public:
  // constructor
  active_face_const_iterator(const Face * p_face,
                             const std::vector<Face> & grid_faces);
  // comparison operator
  bool operator==(const active_face_const_iterator & other) const;
  // comparison operator
  bool operator!=(const active_face_const_iterator & other) const;
  // increment operator
  active_face_const_iterator & operator++();
  // access pointer operator
  inline const Face* operator->() const { return p_face; }
  // access reference operator
  inline const Face& operator*() const { return *p_face; }

 protected:
  void increment_internal_();
  const Face * p_face;
  const std::vector<Face> & m_grid_faces;
};

}  // end namespace mesh
