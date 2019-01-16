#pragma once

#include "Point.hpp"
#include "Plane.hpp"
#include "Polygon.hpp"
#include "Polyhedron.hpp"
#include <Exceptions.hpp>
// #include "PolyGroup.hpp"
#include <CollisionGJK.hpp>
#include <utils.hpp>


/* This module contains various algorithms for simple shape intersections.
 * As opposed to CollisionGJK module, here we actually compute the section
 * data. These algorithms can be combined and utilized to get intersection data
 * of complex shape e.g. a mesh with a fracture.
 *
 * All functions in this module are boolean and return true in the case of
 * collision. The resulting section data is written into the last input argument.
 * If the last aregument is a vector, the data is appended.
 */

namespace angem
{

// get intersection of two planes
// returns true if the planes are not parallel
// the result of the intersection is saved into Line intersection.
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


// collision of a polygon with a plane
// can be 1 points, two points, or zero points
template <typename Scalar>
bool collision(const Polygon<Scalar>        & poly,
               const Plane<Scalar>          & plane,
               std::vector<Point<3,Scalar>> & intersection,
               const double                   tol = 1e-10)
{
  // call collision of all edges
  bool result = false;
  const auto & pts = poly.get_points();
  for (std::size_t i=0; i<pts.size(); ++i)
  {
    bool loc_collision = false;
    if (i < pts.size() - 1)
      loc_collision = collision(pts[i], pts[i+1], plane, intersection, tol);
    else
      loc_collision = collision(pts[i], pts[0], plane, intersection, tol);
    if (loc_collision)
      result = true;
  }

  return result;
}



// collision of a polygon with a plane
// can be 1 points, two points, or zero points
template <typename Scalar>
bool collision(const Polygon<Scalar>        & poly1,
               const Polygon<Scalar>        & poly2,
               std::vector<Point<3,Scalar>> & intersection,
               const double                   tol = 1e-10)
{
  // CollisionGJK<double> collision_gjk;
  // if (!collision_gjk.check(poly1, poly2))
  // {
  //   std::cout << "no GJK" << std::endl;
  //   return false;
  // }
  // std::cout << "yes GJK" << std::endl;

  if (poly1.plane.normal().parallel(poly2.plane.normal(), tol))
  {
    // 1. find vertices of each poly inside another
    // 2. find intersection of edges if any points inside
    PointSet<3,Scalar> pset(tol * 1.5);
    // 1.
    const auto & pts1 = poly1.get_points();
    const auto & pts2 = poly2.get_points();
    bool all_inside1 = true, all_inside2 = true;
    for (const auto & p : pts1)
      if (poly2.point_inside(p, tol))
        pset.insert(p);
      else
        all_inside2 = false;

    for (const auto & p : pts2)
      if (poly1.point_inside(p, tol))
        pset.insert(p);
      else
        all_inside1 = false;

    // std::cout << pset.points << std::endl;
    if (all_inside1)
    {
      pset.points.clear();
      pset.points = pts2;
    }
    if (all_inside2)
    {
      pset.points.clear();
      pset.points = pts1;
    }
    // if (all_inside1 or all_inside2)
    //   throw std::runtime_error("wtf");

    // 2.
    if ( !pset.empty() and !all_inside1 and !all_inside2 )
      for (const auto & edge1 : poly1.get_edges())
      {
        std::vector<Point<3,Scalar>> v_points;
        Plane<Scalar> side = poly1.get_side(edge1);
        for (const auto & edge2 : poly2.get_edges())
        {
          collision(pts2[edge2.first], pts2[edge2.second],
                    side, v_points, tol);
          for (const auto & p : v_points)
            if (poly1.point_inside(p))
              pset.insert(p);
        }
      }

    for (const auto & p: pset.points)
      intersection.push_back(p);

    if (pset.size() > 0)
      return true;
    else
      return false;
  }
  else // two polygons in non-parallel planes
  {
    std::vector<Point<3,Scalar>> v_section;
    angem::collision(poly1, poly2.plane, v_section, tol);
    bool result = false;
    for (const auto & p : v_section)
      if (poly2.point_inside(p, tol) and poly1.point_inside(p, tol))
      {
        intersection.push_back(p);
        result = true;
      }
    return result;
  }


  return true;
}


// intersection of a segment with plane
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


// marks polygons above fracture as 1
// polygons below fracture as 0
template <typename Scalar>
void split(const Polygon<Scalar> & poly,
           const Plane<Scalar>   & plane,
           PolyGroup<Scalar>     & result,
           const int               marker_below = 0,
           const int               marker_above = 1)
{
  std::vector<Point<3,Scalar>> section;
  collision(poly, plane, section);

  if (section.size() == 0)
  {
    std::vector<std::size_t> indices;
    bool above = false;
    for (const auto & p : poly.get_points())
    {
      const std::size_t ind = result.vertices.insert(p);
      indices.push_back(ind);

      if (plane.above(p))  // technically need to check only one
        above = true;
    }

    result.polygons.push_back(indices);

    // assign markers
    if (above)
      result.markers.push_back(marker_above);
    else
      result.markers.push_back(marker_below);

    return;
  }

  std::vector<std::size_t> above, below;

  for (const auto p : poly.get_points())
  {
    const std::size_t ind = result.vertices.insert(p);
    if (plane.above(p))
      above.push_back(ind);
    else
      below.push_back(ind);

  }

  for (Point<3,Scalar> & p : section)
  {
    const std::size_t ind = result.vertices.insert(p);
    above.push_back(ind);
    below.push_back(ind);
  }


  if (above.size() > 2)
  {
    result.polygons.push_back(std::move(above));
    result.markers.push_back(marker_above);
  }
  if (below.size() > 2)
  {
    result.polygons.push_back(std::move(below));
    result.markers.push_back(marker_below);
  }
}


// throws std::runtime_error
template <typename Scalar>
bool collision(const Line<3,Scalar> & line,
               const Plane<Scalar>  & plane,
               Point<3,Scalar>      & intersection)
{
  // Plane : (p - p0) · n = 0
  // line p = d l + l0
  // intersection: d = (p0 - l0) · n / (l · n)
  // Note: if (l · n) == 0 then line is parallel to the plane
  if (line.direction.dot(plane.normal()) < 1e-16)
  {
    if (plane.distance(line.point) < 1e-16)
      throw std::runtime_error("line and plane coinside.");
    return false;
  }

  const Scalar d = (plane.point - line.point).dot(plane.normal()) /
      line.direction.dot(plane.normal());
  intersection = line.point + d*line.direction;
  return true;
}


// section is a vector cause line can reside on polygon
// appends to vector intersection
// note: polygon should have sorted points
template <typename Scalar>
bool collision(const Line<3,Scalar>         & line,
               const Polygon<Scalar>        & poly,
               std::vector<Point<3,Scalar>> & intersection)
{
  // find intersection between polygon plane and line
  Point<3,Scalar> p;
  const bool colinear = collision(line, poly.plane, p);
  if (colinear)
    return false;

  if (poly.point_inside(p), 1e-4)
  {
    intersection.push_back(p);
    return true;
  }
  else
    return false;
}


// collision of a line segment with a polyhedron
template <typename Scalar>
bool collision(const Point<3,Scalar>        & l0,
               const Point<3,Scalar>        & l1,
               const Polyhedron<Scalar>     & poly,
               std::vector<Point<3,Scalar>> & intersection,
               const double                   tol = 1e-10)
{
  std::vector<Point<3,Scalar>> new_section;

  if (poly.point_inside(l0))
    new_section.push_back(l0);
  if (poly.point_inside(l1))
    new_section.push_back(l1);

  if (new_section.size() < 2)
  {
    const auto & points = poly.get_points();
    for (const auto & face : poly.get_faces())
    {
      const std::size_t ibegin = new_section.size();

      Polygon<Scalar> poly_face(points, face);
      collision(l0, l1,  poly_face.plane, new_section, tol);

      for (std::size_t i=ibegin; i<new_section.size(); ++i)
        if (!poly_face.point_inside(new_section[i], tol))
          new_section.erase(new_section.begin() + i);
    }

    remove_duplicates(new_section, tol);
  }

  for (const auto & p : new_section)
    intersection.push_back(p);

  if (new_section.empty())
    return false;
  else
    return true;
}

}  // end namespace
