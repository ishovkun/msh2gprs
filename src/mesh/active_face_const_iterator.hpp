#pragma once

#include "Face.hpp"

namespace mesh {

class active_face_const_iterator
{
 public:
  // constructor
  active_face_const_iterator(const Face * p_face);
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
  const Face * p_face;
};

}  // end namespace mesh
