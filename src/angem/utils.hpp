#pragma once

#include <limits>  // std::numeric_limits
#include <Point.hpp>

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

}
