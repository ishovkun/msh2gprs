#pragma one

#include "Point.hpp"
#include "Plane.hpp"

namespace angem
{

template <typename Scalar>
class Polygon
{
 public:
  Polygon(std::vector< Point<3,Scalar> > & bounding_points);

 protected:
  std::vector< Point<3> > bounding_points;
  Plane plane;
};


template <typename Scalar>
Polygon<Scalar>::Polygon(std::vector< Point<3,Scalar> > & bounding_points)
    :
    bounding_points(bounding_points),
    plane(Point<3,Scalar>(0, 0, 0), Point<3,Scalar>(1, 0, 0), Point<3,Scalar>(0, 1, 0))
{
  assert(bounding_points.size() > 2);

  plane = Plane(bounding_points[0],
                bounding_points[1],
                bounding_points[2]);

}
}
