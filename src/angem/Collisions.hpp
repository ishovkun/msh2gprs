#pragma once

#include "Point.hpp"
#include "Plane.hpp"
#include "Polygon.hpp"
// #include "PolyGroup.hpp"

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
bool collision(const Polygon<Scalar>        & poly,
               const Plane<Scalar>          & plane,
               std::vector<Point<3,Scalar>> & intersection)
{
  // call collision of all edges
  bool result = false;
  const auto & pts = poly.get_points();
  for (std::size_t i=0; i<pts.size(); ++i)
  {
    bool loc_collision = false;
    if (i < pts.size() - 1)
      loc_collision = collision(pts[i], pts[i+1], plane, intersection);
    else
      loc_collision = collision(pts[i], pts[0], plane, intersection);
    if (loc_collision)
      result = true;
  }

  return result;
}


// intersection is appended to!
template <typename Scalar>
bool collision(const Point<3,Scalar>        & l0,
               const Point<3,Scalar>        & l1,
               const Plane<Scalar>          & plane,
               std::vector<Point<3,Scalar>> & intersection,
               const double                   tol = 1e-10)
{
  // Plane : (p - p0) · n = 0
  // line p = d l + l0
  // segment : l0, l1
  // intersection: d = (p0 - l0) · n / (l · n)
  // call collision of all edges
  const Scalar d1 = plane.distance(l0);
  const Scalar d2 = plane.distance(l1);

  // both points are on the plane
  if (fabs(d1) + fabs(d2) < tol)
  {
    intersection.push_back(l0);
    intersection.push_back(l1);
    return true;
  }

  if (d1*d2 > 0)  // both points on one side of plane
    return false;

  // compute intersection point
  const Point<3,Scalar> l = l1 - l0;
  const Scalar d = (plane.point - l0).dot(plane.normal()) /
                    l.dot(plane.normal());
  intersection.push_back(l0 + d * l);
  return true;
}


// template <typename Scalar>
// void split(const Polygon<Scalar> & poly,
//            const Plane<Scalar>   & plane,
//            PolyGroup<Scalar>     & result)
// {
//   std::vector<Point<3,Scalar>> section;
//   collision(poly, plane, section);
//   if (section.size() == 0)
//   {
//     std::vector<Point<3,Scalar>*> p_points;
//     for (const auto p_point : poly.get_points())
//     {
//       result.vertices.push_back(*p_point);
//       auto & back = result.vertices.back();
//       p_points.push_back(&back);
//     }

//     std::cout << "case 1" << std::endl;
//     result.polygons.push_back(Polygon<Scalar>(p_points));
//     return;
//   }

//   std::cout << "case 2" << std::endl;
//   std::vector<Point<3,Scalar>*> points_above, points_below;
//   for (const auto p_point : poly.get_points())
//   {
//     result.vertices.push_back(*p_point);
//     Point<3,Scalar> & p = result.vertices.back();
//     std::cout << "adding " << p << " (" << &p << ")" << std::endl;
//     if (plane.above(p))
//       points_above.push_back(&p);
//     else
//       points_below.push_back(&p);
//   }

//   for (Point<3,Scalar> & p : section)
//   {
//     result.vertices.push_back(p);
//     Point<3,Scalar> & pp = result.vertices.back();
//     std::cout << "adding " << pp << " (" << &pp << ")"<< std::endl;
//     std::cout << "same as " <<  &(result.vertices.back()) << ")"<< std::endl;
//     points_above.push_back(&pp);
//     points_below.push_back(&pp);
//   }


//   // if (points_above.size() > 2)
//   //   result.polygons.push_back(std::move(Polygon<Scalar>(points_above)));
//   // if (points_below.size() > 2)
//   //   result.polygons.push_back(std::move(Polygon<Scalar>(points_below)));

//   // std::cout << "poly_above" << std::endl;
//   // std::cout << result.polygons[0] << std::endl;
//   std::cout << "points_below" << std::endl;
//   for (const auto & pp : points_below)
//     std::cout << *pp << " (" << pp << ")" << std::endl;
//   // std::cout << "poly_below" << std::endl;
//   // std::cout << result.polygons[1] << std::endl;
//   std::cout << "all addresses" << std::endl;
//   for (const auto & p : result.vertices)
//     std::cout << p << " (" << &p << ")"<< std::endl;

// }



}  // end namespace
