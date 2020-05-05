#include "FeValues.hpp"

namespace discretization {

template <>
FeValues<VTK_ID::WedgeID>::FeValues()
    : _center(0.0, 0.0, 0.0)
{}

const double sqrt3 = 1.732050807568877;

template <>
double FeValues<VTK_ID::WedgeID>::eval_(const Point & point, const size_t vertex) const
{
  switch (vertex)  // node numbering in vtk order
  {
    case 0: // φ₁ = ⅙ (1 + 2u) (1 - w)
      return 1.0 / 6.0 * (1 + 2*point[0]) * (1 - point[2]);
    case 2: // φ₂ = ⅙ (1 - u - √3 v) (1 - w)
      return 1.0 / 6.0 * (1 - point[0] - sqrt3*point[1]) * (1 - point[2]);
    case 1:  //φ₃ = ⅙(1 - u + √3 v)(1 - w)
      return 1.0 / 6.0 * (1 - point[0] + sqrt3*point[1]) * (1 - point[2]);
    case 3: // φ₄ = ⅙ (1 + 2u) (1 + w)
      return 1.0 / 6.0 * (1 + 2*point[0]) * (1 + point[2]);
    case 5: // φ₂ = ⅙ (1 - u - √3 v) (1 + w)
      return 1.0 / 6.0 * (1 - point[0] - sqrt3*point[1]) * (1 + point[2]);
    case 4:  //φ₃ = ⅙(1 - u + √3 v)(1 + w)
      return 1.0 / 6.0 * (1 - point[0] + sqrt3*point[1]) * (1 + point[2]);
    default:
      throw std::invalid_argument( "vertex cannot be larger than 5" );
  }
}

template <>
Point FeValues<VTK_ID::WedgeID>::eval_derivative_(const Point & point, const size_t vertex) const
{
  switch (vertex) // node numbering in vtk order
  {
    case 0:
      return {1.0/3.0 * (1 - point[2]),
              0,
              -1.0/6 *(1 + 2*point[0])};
    case 2:
      return {-1.0/6   * (1 - point[2]),
              -sqrt3/6 * (1 - point[2]),
              -1.0/6   * (1 - point[0] - sqrt3*point[1])};
    case 1:
      return {-1.0/6   * (1 - point[2]),
              +sqrt3/6 * (1 - point[2]),
              -1.0/6   * (1 - point[0] + sqrt3*point[1])};
    case 3:
      return {1.0/3 * (1 + point[2]),
              0,
              1.0/6 * (1 + 2*point[0])};
    case 5:
      return {-1.0/6    * (1 + point[2]),
              -sqrt3/6 * (1 + point[2]),
              1.0/6    * (1 - point[0] - sqrt3*point[1])};
    case 4:
      return {-1.0/6  * (1 + point[2]),
              sqrt3/6 * (1 + point[2]),
              1.0/6   * (1 - point[0] + sqrt3*point[1])};
    default:
      throw std::invalid_argument( "vertex cannot be larger than 5" );
  }
}

template <>
std::vector<Point> FeValues<VTK_ID::WedgeID>::get_master_integration_points() const
{
  // 6-point quadrature
  // taken from http://www.softeng.rl.ac.uk/st/projects/felib4/Docs/html/Level-0/qwdg6/qwdg6.html
  return {
    {1, 0, -1},
    {1, 0,  1},
    {-0.5, -0.5*sqrt3, -1},
    {-0.5, -0.5*sqrt3,  1},
    {-0.5,  0.5*sqrt3, -1},
    {-0.5,  0.5*sqrt3,  1},
  };
}

template <>
std::vector<double> FeValues<VTK_ID::WedgeID>::get_master_integration_weights() const
{
  const double w = 0.25 * sqrt3;  // volum=8; 6 points
  return {w, w, w, w, w, w};
}

}  // end namespace discretization
