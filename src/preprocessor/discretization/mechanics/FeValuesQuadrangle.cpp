#include "FeValues.hpp"

namespace discretization {

template <>
FeValues<VTK_ID::QuadrangleID>::FeValues()
    : _center(0.0, 0.0, 0.0)
{}

template <>
double FeValues<VTK_ID::QuadrangleID>::eval_(const Point & x, const size_t vertex) const
{
  switch (vertex)  // φ = ¼ (1 ± ѯ) (1 ± η)
  {
    case 0:
      return 0.25 * (1 - x[0]) * (1 - x[1]);
    case 1:
      return 0.25 * (1 + x[0]) * (1 - x[1]);
    case 2:
      return 0.25 * (1 + x[0]) * (1 + x[1]);
    case 3:
      return 0.25 * (1 - x[0]) * (1 + x[1]);
    default:
      throw std::invalid_argument( "vertex cannot be larger than 3" );
  }
}

template <>
Point FeValues<VTK_ID::QuadrangleID>::eval_derivative_(const Point & x, const size_t vertex) const
{
  switch (vertex)
  {
    case 0:
      return  {-0.25 * (1 - x[1]), -0.25 * (1 - x[0]), 0};
    case 1:
      return  {+0.25 * (1 - x[1]), -0.25 * (1 + x[0]), 0};
    case 2:
      return  {+0.25 * (1 + x[1]), +0.25 * (1 + x[0]), 0};
    case 3:
      return  {-0.25 * (1 + x[1]), +0.25 * (1 - x[0]), 0};
    default:
      throw std::invalid_argument( "vertex cannot be larger than 3" );
  }
}

template <>
std::vector<Point> FeValues<VTK_ID::QuadrangleID>::get_master_integration_points() const
{
  // 3-point rule, stolen from gmsh
  return {
    Point( 0.816496580928,  0.000000000000, 0.000000000000  ),
    Point( -0.408248290464, 0.840896415255, 0.000000000000 ),
    Point( -0.408248290464, -0.840896415255,  0.000000000000 )
  };
}

template <>
std::vector<double> FeValues<VTK_ID::QuadrangleID>::get_master_integration_weights() const
{
  // 3-point rule, stolen from gmsh
  const double w1 = 1.3333333333;
  return {w1, w1, w1};;
}

}  // end namespace discretization
