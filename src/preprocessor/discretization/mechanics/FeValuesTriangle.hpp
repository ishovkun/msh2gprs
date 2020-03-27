#pragma once
#include "FeValues.hpp"
#include <cassert>

namespace discretization {

template<> constexpr size_t N_ELEMENT_VERTICES<VTK_ID::TriangleID> = 3;

template <>
double FeValues<VTK_ID::TriangleID>::eval_(const Point & point, const size_t vertex)
{
  switch (vertex)
  {
    case 0:
      return  1.0 - point[0] - point[1]; // phi_0 = (1 - u - v)
    case 1:
      return point[0];  // phi_1 = u
    case 2:
      return point[1];  // phi_2 = v
    default:
      throw std::invalid_argument( "vertex cannot be larger than 2" );
  }
}

template <>
Point FeValues<VTK_ID::TriangleID>::eval_derivative_(const Point & point, const size_t vertex)
{
  switch (vertex)
  {
    case 0:
      return {-1, -1, 0};
    case 1:
      return {1, 0, 0};
    case 2:
      return {0, 1, 0};
    default:
      throw std::invalid_argument( "vertex cannot be larger than 2" );
  }
}

template <>
std::vector<Point> FeValues<VTK_ID::TriangleID>::get_master_integration_points() const
{
  return {Point(0.125, 0.125, 0)};
}

}  // end namespace discretization
