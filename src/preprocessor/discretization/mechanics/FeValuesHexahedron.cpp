// #pragma once

#include "FeValues.hpp"

namespace discretization {

template <>
FeValues<VTK_ID::HexahedronID>::FeValues()
    : _center(0.0, 0.0, 0.0)
{}

template <>
double FeValues<VTK_ID::HexahedronID>::eval_(const Point & x, const size_t vertex) const
{
  switch (vertex)  //φ = ⅛ (1 ± Ѯ)·(1 ± η)·(1 ± ζ)
  {
    case 0:
      return  0.125 * ( 1 - x[0] ) * ( 1 - x[1] ) * (1 - x[2]);
    case 1:
      return  0.125 * ( 1 + x[0] ) * ( 1 - x[1] ) * (1 - x[2]);
    case 2:
      return  0.125 * ( 1 + x[0] ) * ( 1 + x[1] ) * (1 - x[2]);
    case 3:
      return  0.125 * ( 1 - x[0] ) * ( 1 + x[1] ) * (1 - x[2]);
    case 4:
      return  0.125 * ( 1 - x[0] ) * ( 1 - x[1] ) * (1 + x[2]);
    case 5:
      return  0.125 * ( 1 + x[0] ) * ( 1 - x[1] ) * (1 + x[2]);
    case 6:
      return  0.125 * ( 1 + x[0] ) * ( 1 + x[1] ) * (1 + x[2]);
    case 7:
      return  0.125 * ( 1 - x[0] ) * ( 1 + x[1] ) * (1 + x[2]);
    default:
      throw std::invalid_argument( "vertex cannot be larger than 7" );
  }
}

template <>
Point FeValues<VTK_ID::HexahedronID>::eval_derivative_(const Point & x, const size_t vertex) const
{
  switch (vertex)
  {
    case 0:
      return  {-0.125 * (1 - x[1]) * (1 - x[2]),
               -0.125 * (1 - x[0]) * (1 - x[2]),
               -0.125 * (1 - x[0]) * (1 - x[1])};
    case 1:
      return  {+0.125 * (1 - x[1]) * (1 - x[2]),
               -0.125 * (1 + x[0]) * (1 - x[2]),
               -0.125 * (1 + x[0]) * (1 - x[1])};
    case 2:
      return  {+0.125 * (1 + x[1]) * (1 - x[2]),
               +0.125 * (1 + x[0]) * (1 - x[2]),
               -0.125 * (1 + x[0]) * (1 + x[1])};
    case 3:
      return  {-0.125 * (1 + x[1]) * (1 - x[2]),
               +0.125 * (1 - x[0]) * (1 - x[2]),
               -0.125 * (1 - x[0]) * (1 + x[1])};
    case 4:
      return  {-0.125 * (1 - x[1]) * (1 + x[2]),
               -0.125 * (1 - x[0]) * (1 + x[2]),
               +0.125 * (1 - x[0]) * (1 - x[1])};
    case 5:
      return  {+0.125 * (1 - x[1]) * (1 + x[2]),
               -0.125 * (1 + x[0]) * (1 + x[2]),
               +0.125 * (1 + x[0]) * (1 - x[1])};
    case 6:
      return  {+0.125 * (1 + x[1]) * (1 + x[2]),
               +0.125 * (1 + x[0]) * (1 + x[2]),
               +0.125 * (1 + x[0]) * (1 + x[1])};
    case 7:
      return  {-0.125 * (1 + x[1]) * (1 + x[2]),
               +0.125 * (1 - x[0]) * (1 + x[2]),
               +0.125 * (1 - x[0]) * (1 + x[1])};
    default:
      throw std::invalid_argument( "vertex cannot be larger than 7" );
  }
}

template <>
std::vector<Point> FeValues<VTK_ID::HexahedronID>::get_master_integration_points() const
{
  // 6-point quadrature
  // Values taken from GMsh/numeric/GaussQuadratureHex.cpp
  return {
    Point(  0.40824826,  0.70710678, -0.57735027 ),
    Point(  0.40824826, -0.70710678, -0.57735027 ),
    Point( -0.40824826,  0.70710678,  0.57735027 ),
    Point( -0.40824826, -0.70710678,  0.57735027 ),
    Point( -0.81649658,  0.00000000, -0.57735027 ),
    Point(  0.81649658,  0.00000000,  0.57735027 )
  };
}

template <>
std::vector<double> FeValues<VTK_ID::HexahedronID>::get_master_integration_weights() const
{
  const double w1 = 8.0 / 6.0;  // volum=8; 6 points
  return {w1, w1, w1, w1, w1, w1};;
}

}  // end namespace discretization
