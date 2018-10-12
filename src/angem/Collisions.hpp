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

  // check that intersection point is within the polygon
  // algorithm: if section point is on the same side of the faces as the
  // mass center, then the point is inside of the polygon
  const auto & poly_verts = poly.get_points();
  Point<3,Scalar> cm = compute_center_mass(poly_verts);
  const auto & normal = poly.plane.normal();

  bool inside = true;
  for (std::size_t i=0; i<poly_verts.size(); ++i)
  {
    const Point<3,Scalar> & v1 = poly_verts[i];
    Point<3,Scalar> v2;
    if (i < poly_verts.size() - 1)
      v2 = poly_verts[i+1];
    else
      v2 = poly_verts[0];

    Point<3,Scalar> p_perp = v1 + normal * (v2 - v1).norm();
    Plane<Scalar> side(v1, v2, p_perp);
    if (side.above(p) != side.above(cm))
      return false;
  }

  intersection.push_back(p);
  return true;
}

}  // end namespace
