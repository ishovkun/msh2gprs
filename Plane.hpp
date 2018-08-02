// Geometrical plane object in 3D

#pragma once

#include "Point.hpp"

namespace angem
{

template <typename Scalar>
class Plane
{
 public:
  // create plane from a point on the plane and a normal vector
  Plane(const Point<3> & point,
        const Point<3> & normal);
  // strike and dip in degrees
  Plane(const Point<3> & point,
        const double     dip_angle,      // -90 <= dip <= 90
        const double     strike_angle);
  // create plane from 3 points
  Plane(const Point<3>   & p1,
        const Point<3>   & p2,
        const Point<3>   & p3);
  friend Point<3> project(const Point<3> & p,
                          const Plane    & plane);

  friend double distance(const Point<3> & p,
                         const Plane    & plane);

  friend bool above(const Point<3> & p,
                    const Plane    & plane);

 protected:
  Point<3> point, normal;
};


}  // end namespace
