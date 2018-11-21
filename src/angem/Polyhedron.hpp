#pragma once

#include <Shape.hpp>
#include <Plane.hpp>
#include <PointSet.hpp>
#include <typeinfo>

namespace angem
{

template<typename Scalar>
class Polyhedron: public Shape<Scalar>
{
 public:
  Polyhedron(const int vtk_id = -1);
  Polyhedron(const std::vector<Point<3,Scalar>>          & vertices,
             const std::vector<std::vector<std::size_t>> & faces,
             const int                                     vtk_id = -1);

  void set_data(const std::vector<Point<3,Scalar>>          & vertices,
                const std::vector<std::vector<std::size_t>> & faces);
  int id() const {return vtk_id;}

  const std::vector<std::vector<std::size_t>> & get_faces() const;
  std::vector<std::vector<std::size_t>> & get_faces();


 protected:
  std::vector<std::vector<std::size_t>> faces;
  const int vtk_id;
};


template<typename Scalar>
Polyhedron<Scalar>::Polyhedron(const int vtk_id)
    :
    vtk_id(vtk_id)
{}


template<typename Scalar>
Polyhedron<Scalar>::Polyhedron(const std::vector<Point<3,Scalar>>          & vertices,
                               const std::vector<std::vector<std::size_t>> & faces,
                               const int vtk_id)
    :
    vtk_id(vtk_id)
{
  set_data(vertices, faces);
}


template<typename Scalar>
void
Polyhedron<Scalar>::set_data(const std::vector<Point<3,Scalar>>          & vertices,
                             const std::vector<std::vector<std::size_t>> & faces)
{
  assert(vertices.size() > 3);
  PointSet<3,Scalar> pset(distance(vertices[0], vertices[1]));
  std::size_t iface = 0;
  for (const auto & face : faces)
  {
    this->faces[iface].reserve(face.size());
    for(const auto ivert : face)
    {
      const Point<3,Scalar> p = vertices[ivert];
      this->faces[iface].push_back(pset.insert(p));
    }
    iface++;
  }

  this->points.reserve(vertices.size());
  for (const auto & p : pset.points)
    this->points.push_back(p);
}


template<typename Scalar>
std::vector<std::vector<std::size_t>> &
Polyhedron<Scalar>::get_faces()
{
  return faces;
}


template<typename Scalar>
const std::vector<std::vector<std::size_t>> &
Polyhedron<Scalar>::get_faces() const
{
  return faces;
}

}
