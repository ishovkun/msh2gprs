// Geometrical plane object in 3D

#pragma once

#include "Point.hpp"

namespace angem
{

class Plane
{
 public:
  Plane(Point<3> & point,
        Point<3> & normal);

  // strike and dip in degrees
  Plane(Point<3>   & point,
        const double dip_angle,      // -90 <= dip <= 90
        const double strike_angle);

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
