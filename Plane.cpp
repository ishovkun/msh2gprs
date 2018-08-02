#include "Plane.hpp"

#include <math.h>    // M_PI

namespace angem
{

Plane::Plane(const Point<3> & point,
             const Point<3> & normal)
    :
    point(point),
    normal(normal)
{
  assert(fabs(normal.norm() - 1) < 1e-12);
}


Plane::Plane(const Point<3> & point,
             const double     dip_angle,
             const double     strike_angle)
    :
    point(point)
{
  assert(dip_angle >= - 90 and dip_angle <= 90);

  double rdip    = dip_angle * M_PI / 180.;
  double rstrike = strike_angle * M_PI / 180.;

  // dip = polar angle between normal and z
  normal[0] = sin(rdip) * cos(rstrike + M_PI/2.);
  normal[1] = sin(rdip) * sin(rstrike + M_PI/2.);
  normal[2] = cos(rdip);
}


Plane::Plane(const Point<3>   & p1,
             const Point<3>   & p2,
             const Point<3>   & p3)
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


double distance(const Point<3> & p,
                const Plane    & plane)
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


bool above(const Point<3> & p,
           const Plane    & plane)
{
  if (distance(p, plane) > 0)
    return true;
  else
    return false;
}


}  // end namespace
