#pragma once

#include "Point.hpp"
#include "Plane.hpp"
#include "Polygon.hpp"
#include "utils.hpp"

namespace angem
{

template <typename Scalar>
class Cell
{
 public:
  Cell(const std::vector<Point<3,Scalar>> & vertices);
  bool intersects(const Polygon<Scalar> poly);


  std::vector<Point<3,Scalar>> vertices;
  Point<3, Scalar> center;
};


template <typename Scalar>
Cell<Scalar>::Cell(const std::vector<Point<3,Scalar>> & vertices)
    :
    vertices(vertices)
{
  assert(vertices.size() > 3);
  assert(!align_on_plane(vertices));
  center = compute_center_mass(vertices);
}


template <typename Scalar>
bool Cell<Scalar>::intersects(const Polygon<Scalar> poly)
{
  /* Algorithm:
   * check if distance between cell center and
   */
  return true;
}

}  // end namespace
