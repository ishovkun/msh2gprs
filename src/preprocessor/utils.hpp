#include <limits>  // std::numeric_limits
#include "Point.hpp"

namespace angem
{

template<int dim, typename Scalar>
Point<dim,Scalar>
compute_center_mass(const std::vector<Point<dim,Scalar>> & points)
{
  Point<dim, Scalar> center = {0, 0, 0};
  for (const auto & p : points)
    center += p;
  center /= static_cast<Scalar>(points.size());
  return center;
}


template<int dim, typename Scalar>
Point<dim,Scalar>
fathest_point_in_direction(const std::vector<Point<dim,Scalar>> & pts,
                           const Point<dim,Scalar> & dir)
{
  Scalar max_dist = - std::numeric_limits<Scalar>::max();
  std::size_t ind = pts.size();

  for (std::size_t i=0; i<pts.size(); ++i)
  {
    Scalar dist = pts[i].dot(dir);
    if (dist > max_dist)
    {
      max_dist = dist;
      ind = i;
    }
  }

  assert( ind != pts.size() );
  return pts[ind];
}


}
