#include "active_face_const_iterator.hpp"

namespace mesh {

active_face_const_iterator::active_face_const_iterator(const Face * p_face)
    : p_face(p_face)
{
  // if (p_face)
    // assert (p_face->is_active() &&
    //         "Trying to get active face iterator for inactive face");
}

}  // end namespace mesh
