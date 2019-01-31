#include <face_iterator.hpp>
#include <mesh_methods.hpp>

namespace mesh
{

face_iterator::
face_iterator(const FaceMap::iterator            & it,
              angem::PointSet<3,double>          & vertices)
    :
    face_it(it),
    mesh_vertices(vertices)
{}


face_iterator::
face_iterator(const face_iterator & other)
    :
    face_it(other.face_it),
    mesh_vertices(other.mesh_vertices)
{}


face_iterator &
face_iterator::operator=(const face_iterator & other)
{
  face_it = other.face_it;
  mesh_vertices = other.mesh_vertices;
  return (*this);
}


bool face_iterator::operator==(const face_iterator & other) const
{
  if (face_it != other.face_it)
    return false;

  return true;
}


bool face_iterator::operator!=(const face_iterator & other) const
{
  return !(*this == other);
}


face_iterator & face_iterator::operator++()
{
  face_it++;
  return (*this);
}


int face_iterator::marker() const
{
  return face_it->second.marker;
}


std::vector<Point> face_iterator::vertices() const
{
  std::vector<std::size_t> ivertices(invert_hash(face_it->first));
  return get_vertex_coordinates(mesh_vertices, ivertices);
}


std::vector<std::size_t> face_iterator::vertex_indices() const
{
  return face_it->second.ordered_indices;
}



}
