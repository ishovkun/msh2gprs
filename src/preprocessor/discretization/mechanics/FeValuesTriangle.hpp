#pragma once
#include "FeValues.hpp"
#include <cassert>

namespace discretization {

template<> constexpr size_t N_ELEMENT_VERTICES<VTK_ID::TriangleID> = 3;

template <>
FeValues<VTK_ID::TriangleID>::FeValues(const mesh::Mesh & grid)
    : _grid(grid)
{}

template <>
double FeValues<VTK_ID::TriangleID>::eval_(const Point & point, const size_t vertex)
{
  assert( vertex < 3 );
  switch (vertex)
  {
    case 0:
      return  1.0 - point[0] - point[1]; // phi_0 = (1 - x - y)
    case 1:
      return point[0];  // phi_1 = x
    case 2:
      return point[1];  // phi_2 = y
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
