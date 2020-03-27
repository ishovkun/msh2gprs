#pragma once
#include "FeValues.hpp"
#include <cassert>

namespace discretization {

template<> constexpr size_t N_ELEMENT_VERTICES<VTK_ID::TetrahedronID> = 4;

template <>
double FeValues<VTK_ID::TetrahedronID>::eval_(const Point & point, const size_t vertex) const
{
  switch (vertex)
  {
    case 0:
      return  1.0 - point[0] - point[1] - point[2]; // phi_0 = (1 - u - v - w)
    case 1:
      return point[0];  // phi_1 = u
    case 2:
      return point[1];  // phi_2 = v
    case 3:
      return point[2];  // phi_2 = w
    default:
      throw std::invalid_argument( "vertex cannot be larger than 3" );
  }
}

template <>
Point FeValues<VTK_ID::TetrahedronID>::eval_derivative_(const Point & point, const size_t vertex) const
{
  switch (vertex)
  {
    case 0:
      return {-1, -1, -1};
    case 1:
      return {1, 0, 0};
    case 2:
      return {0, 1, 0};
    case 3:
      return {0, 0, 1};
    default:
      throw std::invalid_argument( "vertex cannot be larger than 3" );
  }
}

template <>
std::vector<Point> FeValues<VTK_ID::TetrahedronID>::get_master_integration_points() const
{
  // return {Point(0.125, 0.125, 0.125)};
  return {Point(0.25, 0.25, 0.25)};
}


}  // end namespace discretization
