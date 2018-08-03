#pragma once

#include "Point.hpp"
#include "Plane.hpp"
#include "Polygon.hpp"

namespace angem
{

template <typename Scalar>
bool collision(const Plane<Scalar> & pl1,
               const Plane<Scalar> & pl2,
               Line<3,Scalar>      & intersection)
{
  /* The method selects a third plane P3 with
   * an implicit equation n3 · P = 0
   where n3 = n1 x n2 and d3 = 0 (meaning it passes through the origin).
   This always works since: (1) L is perpendicular to P3 and thus
   intersects it, and (2) the vectors n1, n2, and n3 are linearly independent.
   Thus the planes P1, P2 and P3 intersect in a unique point P0 which must be on L.
   Using the formula for the intersection of 3 planes,
   where d3 = 0 for P3, we get:
        (d2 n1 - d1 n2) x n3
   p3 = --------------------
           (n1 x n2) · n3
   Ref:
   http://geomalgorithms.com/a05-_intersect-1.html
  */

  // direction of intersection line n 3 = n1 x n2
  Point<3,Scalar> n3 = pl1.normal().cross(pl2.normal());
  if (n3.norm() < 1e-16)
    return false;
  n3.normalize();


  intersection.direction = n3;

  const auto & n1 = pl1.normal();
  const auto & d1 = pl1.d;

  const auto & n2 = pl2.normal();
  const auto & d2 = pl2.d;

  Point<3,Scalar> numerator = (d2*n1 - d1*n2).cross(n3);
  Scalar denumenator = (n1.cross(n2)).dot(n3);

  intersection.point = numerator / denumenator;

  return true;
}


template <typename Scalar>
bool collision(const Polygon<2,Scalar> & poly,
               const Line<2,Scalar>    & line)
{

}


template <typename Scalar>
bool collision(const Polygon<3,Scalar> & poly1,
               const Polygon<3,Scalar> & poly2)
{
  /*
   * Algorithm:
   * 1. find intersection of planes and exit if parallel
   * 2. find out if intersects line lies within each polygon
   *    a. Make Polygon1 and intersection plane 2D objects
   *    b. Determine if they collide in 2D
   */

  // 1.
  Line<3,Scalar> line_section;
  const bool non_parallel = collision(poly1.plane, poly2.plane, line_section);
  if (!non_parallel)
    return false;

  // 2.
  // find local coordinates of intersection line
  const Point<3,Scalar> lp = poly1.plane.local_coordinates(line_section.point);
  const Point<3,Scalar> ldirection = poly1.plane.local_coordinates(line_section.point);
  // cast to 2d
  Line<2,Scalar> section2d;
  section2d.point = Point<2,Scalar>(lp(0), lp(1));
  section2d.direction = Point<2,Scalar>(ldirection(0), ldirection(1));

  // find local coordinates of bounding points
  Polygon<2,Scalar> poly2d = poly1.to_2D();

  return collision(poly2d, section2d);
}

}  // end namespace
