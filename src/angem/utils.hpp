#pragma once

#include <limits>  // std::numeric_limits
#include <Point.hpp>
#include <PointSet.hpp>

namespace angem
{

template<int dim, typename Scalar>
Point<dim,Scalar>
compute_center_mass(const std::vector<Point<dim,Scalar>> & points)
{
  Point<dim, Scalar> center = {0, 0, 0};
  for (const auto & p : points)
    center += p;
  center /= static_cast<Scalar>(points.size());
  return center;
}


template<int dim, typename Scalar>
Point<dim,Scalar>
compute_center_mass(const std::vector<Point<dim,Scalar> *> & points)
{
  Point<dim, Scalar> center = {0, 0, 0};
  for (const auto & p : points)
    center += *p;
  center /= static_cast<Scalar>(points.size());
  return center;
}


// this function is O(n²)
template<int dim, typename Scalar>
void
remove_duplicates(const std::vector<Point<dim,Scalar>> & points,
                  std::vector<Point<dim,Scalar>>       & result,
                  const double                           tolerance = 0)
{
  // returns result vector that contains only unique entries of points vector
  // two points are considered duplicate if the distance between them is
  // less than tolerance

  if (result.size() != 0)
    result.clear();

  for (const auto & p : points)
    if (find(p, result, tolerance) == result.size())
      result.push_back(p);
}


// remove duplicates (with tolerance) from the vector of points
// this function is O(n)
template<int dim, typename Scalar>
void remove_duplicates(std::vector<Point<dim,Scalar>> & points,
                       const double tolerance = 1e-6)
{
  PointSet<dim,Scalar> pset;
  for (const auto & p : points)
    pset.insert(p);
  points = std::move(pset.points);
}


template<int dim, typename Scalar>
Scalar triangle_area(const Point<dim,Scalar> & p1,
                     const Point<dim,Scalar> & p2,
                     const Point<dim,Scalar> & p3)
{
  return static_cast<Scalar>(0.5) * (cross( p2 - p1, p3 - p1 )).norm();
}


// convert radians to degrees
template<typename Scalar>
inline
Scalar degrees(const Scalar angle)
{
  return angle / static_cast<Scalar>(M_PI) * static_cast<Scalar>(180.);
}


// convert degrees to degrees
template<typename Scalar>
inline
Scalar radians(const Scalar angle)
{
  return angle * static_cast<Scalar>(M_PI) / static_cast<Scalar>(180.);
}

}
