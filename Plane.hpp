// Geometrical plane object in 3D
#pragma once

#include "Point.hpp"

#include <math.h>    // M_PI

namespace angem
{

template <typename Scalar>
class Plane
{
 public:
  // create plane from a point on the plane and a normal vector
  Plane(const Point<3,Scalar> & point,
        const Point<3,Scalar> & normal);
  // strike and dip in degrees
  Plane(const Point<3,Scalar> & point,
        const Scalar          & dip_angle,      // -90 <= dip <= 90
        const Scalar          & strike_angle);
  // create plane from 3 points
  Plane(const Point<3,Scalar> & p1,
        const Point<3,Scalar> & p2,
        const Point<3,Scalar> & p3);

  template <typename T>
  friend Point<3,T> project(const Point<3,T> & p,
                            const Plane<T>   & plane);

  template <typename T>
  friend T distance(const Point<3,T> & p,
                    const Plane<T>   & plane);

  template <typename T>
  friend bool above(const Point<3,T> & p,
                    const Plane<T>   & plane);

 protected:
  Point<3,Scalar> point, normal;
};


template <typename Scalar>
Plane<Scalar>::Plane(const Point<3,Scalar> & point,
                     const Point<3,Scalar> & normal)
    :
    point(point),
    normal(normal)
{
  assert(fabs(normal.norm() - 1) < 1e-12);
}


template <typename Scalar>
Plane<Scalar>::Plane(const Point<3,Scalar> & point,
                     const Scalar          & dip_angle,
                     const Scalar          & strike_angle)
    :
    point(point)
{
  assert(dip_angle >= - 90 and dip_angle <= 90);

  Scalar rdip    = dip_angle * M_PI / 180.;
  Scalar rstrike = strike_angle * M_PI / 180.;

  // dip = polar angle between normal and z
  normal[0] = sin(rdip) * cos(rstrike + M_PI/2.);
  normal[1] = sin(rdip) * sin(rstrike + M_PI/2.);
  normal[2] = cos(rdip);
}


template <typename Scalar>
Plane<Scalar>::Plane(const Point<3,Scalar> & p1,
                     const Point<3,Scalar> & p2,
                     const Point<3,Scalar> & p3)
    :
    point(p1)
{
  assert(p1 != p2 and p2 != p3 and p1 != p3);
  // define two tangent vectors
  const Point<3> t1 = p1 - p2;
  const Point<3> t2 = p1 - p3;
  normal = cross_product(t1, t2);
  normal.normalize();
}


template <typename Scalar>
Scalar distance(const Point<3,Scalar> & p,
                const Plane<Scalar>   & plane)
{
  /* dot product of point by perpendicular vector:
   * if < 0: point is below the plane
   */
  /* signed distance from point p (vertex) to plane
   * with normal n and containing point x0:
   * d = dot(p - x0, n)
   * if d<0, point below the plane
   */
  double result = 0;
  for (std::size_t i=0; i<3; ++i)
    result += (p(i) - plane.point(i)) * plane.normal(i);
  return result;
}


template <typename Scalar>
bool above(const Point<3,Scalar> & p,
           const Plane<Scalar>   & plane)
{
  if (distance(p, plane) > 0)
    return true;
  else
    return false;
}


template <typename Scalar>
bool align_on_plane(const std::vector<Point<3, Scalar>> & points)
{
  assert(points.size() > 2);

  if (points.size() == 3)
    return true;

  Plane<Scalar> plane(points[0], points[1], points[2]);
  for (std::size_t i=3; i<points.size(); ++i)
    if ( fabs(distance(points[i], plane)) > 1e-16 )
      return false;

  return true;
}


}  // end namespace
