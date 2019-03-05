#pragma once

#include <Shape.hpp>
#include <Plane.hpp>
#include <PointSet.hpp>
#include <Polygon.hpp>
#include <PolyGroup.hpp>
#include <typeinfo>
#include <exception>

namespace angem
{

template<typename Scalar>
class Polyhedron: public Shape<Scalar>
{
 public:
  // Constructors
  Polyhedron(const int vtk_id = -1);
  Polyhedron(const std::vector<Point<3,Scalar>>          & vertices,
             const std::vector<std::vector<std::size_t>> & faces,
             const int                                     vtk_id = -1);
  // setter
  void set_data(const std::vector<Point<3,Scalar>>          & vertices,
                const std::vector<std::vector<std::size_t>> & faces);
  // getters
  int id() const {return vtk_id;}
  virtual Scalar volume() const;
  bool point_inside(const Point<3,Scalar> & p) const;

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
  this->faces.resize(faces.size());

  PointSet<3,Scalar> pset;
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


template<typename Scalar>
Scalar Polyhedron<Scalar>::volume() const
{
  const Point<3,Scalar> c = this->center();
  Scalar vol = 0;
  for (const auto & face_indices : get_faces())
  {
    const auto face_poly = Polygon<Scalar>(this->points, face_indices);
    const Scalar face_area = face_poly.area();
    const Scalar h = c.distance(face_poly.plane.project_point(c));
    vol += 1./3. * h * face_area;
  }
  return vol;
}

// template<typename Scalar>
// Scalar Polyhedron<Scalar>::volume() const
// {
//   // std::cout << "wrong class" << std::endl;
//   // abort();
//   // i don't know how it works
//   // ask mohammad
//   const std::size_t n_faces = faces.size();
//   std::vector<Scalar> face_area(n_faces);
//   // no idea what those are
//   std::vector<Point<3,Scalar>> FG(n_faces), Fn(n_faces);

//   // compute some face area quantities
//   for (std::size_t i=0; i<n_faces; ++i)
//   {
//     const auto & face = faces[i];
//     const Point<3,Scalar> p0 = this->points[face[0]];
//     for (std::size_t j=1; j<face.size(); ++j)
//     {
//       const Point<3,Scalar> & p1 = this->points[face[j]];
//       const Point<3,Scalar> & p2 = this->points[face[j+1]];
//       const Point<3,Scalar> u = p1 - p0;
//       const Point<3,Scalar> v = p2 - p0;
//       const Scalar nx = (u[0]*v[2] - v[1]*u[0]);
//       const Scalar ny = (v[0]*u[2] - u[0]*v[2]);
//       const Scalar nz = (u[0]*v[1] - u[1]*v[0]);
//       const Scalar area_tmp = .5*sqrt(nx*nx + ny*ny + nz*nz);

//       face_area[i] += area_tmp;

//       FG[i] += area_tmp * (p0 + p1 + p2) / 3.;
//       Fn[i][0] += 0.5 * nx;
//       Fn[i][1] += 0.5 * ny;
//       Fn[i][2] += 0.5 * nz;
//     }

//     FG[i] /= face_area[i];
//     Fn[i] /= face_area[i];

//     Scalar nl = std::sqrt(Fn[i].norm());
//     Fn[i] /= nl;
//   }

//   // finally compute volume
//   Scalar vol = 0;
//   const Point<3,Scalar> c = this->center();
//   for (std::size_t j=0; j<n_faces; ++j)
//   {
//     Scalar h;
//     for (int d = 0; d < 3; ++d)
//       h += Fn[j][d]*(FG[j][d] - c[d]);

//     vol += fabs(h*face_area[j]) / 3.;
//   }

//   return vol;
// }


template<typename Scalar>
bool Polyhedron<Scalar>::point_inside(const Point<3,Scalar> & p) const
{
  const Point<3,Scalar> c = this->center();
  for (const auto face : faces)
  {
    Plane<Scalar> plane(this->points[face[0]],
                        this->points[face[1]],
                        this->points[face[2]]);
    if (plane.above(p) != plane.above(c))
      return false;
  }
  return true;
}


template<typename Scalar>
std::ostream &operator<<(std::ostream             & os,
                         const Polyhedron<Scalar> & poly)
{
  const auto & points = poly.get_points();
  const auto & faces = poly.get_faces();
  os << points.size() << " vertices:" << std::endl;
  os << points;
  os << faces.size() << " faces:" << std::endl;
  for (const auto & face : faces)
  {
    for (const auto ivertex: face)
      os << ivertex << "\t";
    os << std::endl;
  }
  return os;
}


}  // end namespace
