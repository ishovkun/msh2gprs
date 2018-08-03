#pragma once

#include "Point.hpp"
#include "Plane.hpp"
#include "Polygon.hpp"
#include "utils.hpp"

namespace angem
{

// template <typename Scalar>
// class Polyhedron
// {
//  public:
//   Polyhedron(const std::vector<Point<3,Scalar>> & vertices);
//   bool intersects(const Polygon<Scalar> poly);
//   const Point<3,Scalar> & center() {return mass_center};


//  private:
//   std::vector<Point<3,Scalar> * > vertices;
//   std::vector<std::vector<std::size_t> > faces;
//   Point<3, Scalar> mass_center;
// };


// template <typename Scalar>
// Polyhedron<Scalar>::Polyhedron(const std::vector<Point<3,Scalar>> & vertices)
//     :
//     vertices(vertices)
// {
//   assert(vertices.size() > 3);
//   assert(!align_on_plane(vertices));
//   center = compute_center_mass(vertices);
// }


// template <typename Scalar>
// bool Polyhedron<Scalar>::intersects(const Polygon<Scalar> poly)
// {
//   /* Algorithm:
//    * check if distance between cell center and
//    */
//   return true;
// }

}  // end namespace
