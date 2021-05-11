#include "FeValues.hpp"
#include <cassert>

namespace discretization {

template <>
FeValues<VTK_ID::TriangleID>::FeValues()
    : _center(0.125, 0.125, 0)
{}

template <>
double FeValues<VTK_ID::TriangleID>::eval_(const Point & point, const size_t vertex) const
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
Point FeValues<VTK_ID::TriangleID>::eval_derivative_(const Point & point, const size_t vertex) const
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
  return {Point(1.0/3.0, 1.0/3.0, 0)};
}

template <>
std::vector<double> FeValues<VTK_ID::TriangleID>::get_master_integration_weights() const
{
  return {0.5};
}

template <>
Point FeValues<VTK_ID::TriangleID>::map_real_to_local_(const Point & xyz) const
{
  /**
   * This code is copied from gmsh and sligtly modified to adapt to the
   * angem Point class.
   */
  const Point & O = _vertex_coord[0];
  const Point d = xyz - O;
  const Point d1 = _vertex_coord[1] - O;
  const Point d2 = _vertex_coord[2] - O;

  const double Jxy = d1[0] * d2[1] - d1[1] * d2[0];
  const double Jxz = d1[0] * d2[2] - d1[2] * d2[0];
  const double Jyz = d1[1] * d2[2] - d1[2] * d2[1];

  Point uvw;

  if((std::abs(Jxy) > std::abs(Jxz)) && (std::abs(Jxy) > std::abs(Jyz))) {
    uvw[0] = (d[0] * d2[1] - d[1] * d2[0]) / Jxy;
    uvw[1] = (d[1] * d1[0] - d[0] * d1[1]) / Jxy;
  }
  else if(std::abs(Jxz) > std::abs(Jyz)) {
    uvw[0] = (d[0] * d2[2] - d[2] * d2[0]) / Jxz;
    uvw[1] = (d[2] * d1[0] - d[0] * d1[2]) / Jxz;
  }
  else {
    uvw[0] = (d[1] * d2[2] - d[2] * d2[1]) / Jyz;
    uvw[1] = (d[2] * d1[1] - d[1] * d1[2]) / Jyz;
  }
  uvw[2] = 0.0;

  return uvw;
}


}  // end namespace discretization