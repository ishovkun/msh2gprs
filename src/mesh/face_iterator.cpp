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
  return std::move(get_vertex_coordinates(mesh_vertices, ivertices));
}


std::vector<std::size_t> face_iterator::vertex_indices() const
{
  // return std::move(invert_hash(face_it->first));
  return face_it->second.ordered_indices;
}



}
