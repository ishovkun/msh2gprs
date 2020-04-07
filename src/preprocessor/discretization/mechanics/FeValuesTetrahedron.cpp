// #pragma once
#include "FeValues.hpp"
#include <cassert>

namespace discretization {

template <>
double FeValues<VTK_ID::TetrahedronID>::eval_(const Point & point, const size_t vertex) const
{
  switch (vertex)
  {
    case 0:
      return  1.0 - point[0] - point[1] - point[2]; // phi_0 = (1 - u - v - w)
    case 1:
      return point[0];                              // phi_1 = u
    case 2:
      return point[1];                              // phi_2 = v
    case 3:
      return point[2];                              // phi_2 = w
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
  return {Point(0.25, 0.25, 0.25)};
}

template <>
std::vector<double> FeValues<VTK_ID::TetrahedronID>::get_master_integration_weights() const
{
  return {1.0};
}

template <>
Point FeValues<VTK_ID::TetrahedronID>::map_real_to_local_(const Point & xyz) const
{
  const angem::Tensor2<3,double> mat = {
    _vertex_coord[1].x() - _vertex_coord[0].x(),
    _vertex_coord[2].x() - _vertex_coord[0].x(),
    _vertex_coord[3].x() - _vertex_coord[0].x(),
    _vertex_coord[1].y() - _vertex_coord[0].y(),
    _vertex_coord[2].y() - _vertex_coord[0].y(),
    _vertex_coord[3].y() - _vertex_coord[0].y(),
    _vertex_coord[1].z() - _vertex_coord[0].z(),
    _vertex_coord[2].z() - _vertex_coord[0].z(),
    _vertex_coord[3].z() - _vertex_coord[0].z()
  };
  const angem::Tensor2<3,double> mat_inv = invert(mat);
  return mat_inv * xyz;
}

}  // end namespace discretization
