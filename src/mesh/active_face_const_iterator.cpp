#include "active_face_const_iterator.hpp"

namespace mesh {

active_face_const_iterator::active_face_const_iterator(const Face * p_face,
                                                       const std::vector<Face> & grid_faces)
    : p_face(p_face), m_grid_faces(grid_faces)
{
  if (p_face)
    assert (p_face->is_active() &&
            "Trying to get active face iterator for inactive face");
}

bool active_face_const_iterator::operator==(const active_face_const_iterator & other) const
{
  if (p_face && other.p_face)
    return *p_face == *other.p_face;
  else
    return p_face == other.p_face;
}

bool active_face_const_iterator::operator!=(const active_face_const_iterator & other) const
{
  return !operator==(other);
}

active_face_const_iterator &
active_face_const_iterator::operator++()
{
  if (p_face->index() == m_grid_faces.size() - 1)
  {
    p_face = nullptr;
    return *this;
  }

  increment_internal_();
  while ( !p_face->is_active() )
  {
    if (p_face->index() == m_grid_faces.size() - 1)
    {
      p_face = nullptr;
      return *this;
    }
    increment_internal_();
  }
  return *this;
}

void active_face_const_iterator::increment_internal_()
{
  const size_t index_next = p_face->index() + 1;
  p_face = &m_grid_faces[index_next];
}

}  // end namespace mesh
