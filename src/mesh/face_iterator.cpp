#include <face_iterator.hpp>
#include <mesh_methods.hpp>

namespace mesh
{

face_iterator::
face_iterator(const FaceMap::iterator            & it,
              angem::PointSet<3,double>          & vertices)
    :
    face_it(it),
    p_mesh_vertices(&vertices)
{}


face_iterator::
face_iterator(const face_iterator & other)
    :
    face_it(other.face_it),
    p_mesh_vertices(other.p_mesh_vertices)
{}


face_iterator &
face_iterator::operator=(const face_iterator & other)
{
  face_it = other.face_it;
  p_mesh_vertices = other.p_mesh_vertices;
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
  return get_vertex_coordinates(p_mesh_vertices, ivertices);
}


std::vector<std::size_t> face_iterator::vertex_indices() const
{
  return face_it->second.ordered_indices;
}


Point face_iterator::normal() const
{
  const auto poly = angem::Polygon<double>(vertices());
  return poly.plane.normal();
}


angem::Polygon<double> face_iterator::polygon() const
{
  return angem::Polygon<double>(vertices());
}


angem::Point<3,double> face_iterator::center() const
{
  const auto verts = vertex_indices();
  angem::Point<3,double> c = {0, 0, 0};
  for (const size_t vert : verts)
    c += (*p_mesh_vertices)[vert];

  for (int comp = 0; comp < 3; comp++)
    c[comp] /= verts.size();
  return c;
}

}
