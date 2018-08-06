#pragma once

#include "Point.hpp"

namespace angem
{

template<typename Scalar>
class Shape
  {
   public:
    // constructors
    Shape();
    Shape(std::vector<Point<3,Scalar>> & point_list);
    // setter
    void set_data(std::vector<Point<3,Scalar>> & point_list);
    // getter
    std::vector<Point<3,Scalar> *> & get_points();
    // check if empty
    bool empty() const;
    // support function for gjk collision algorithm
    virtual Point<3,Scalar> support(Point<3, Scalar> & direction) const;

   protected:
    std::vector<Point<3,Scalar> *> points;
};


template<typename Scalar>
Shape<Scalar>::Shape()
{}


template<typename Scalar>
Shape<Scalar>::Shape(std::vector<Point<3,Scalar>> & point_list)
{
  set_data(point_list);
}


template<typename Scalar>
void
Shape<Scalar>::set_data(std::vector<Point<3,Scalar>> & point_list)
{
  points.reserve(point_list.size());
  for (auto & p : point_list)
    points.push_back(&p);
}


template<typename Scalar>
bool
Shape<Scalar>::empty() const
{
  if (points.size() > 0)
    return false;
  else
    return true;
}


template<typename Scalar>
Point<3,Scalar>
Shape<Scalar>::support(Point<3,Scalar> & direction) const
{
  assert(!empty());

  // get the farthest point in direction p
  Scalar max_dist = - std::numeric_limits<Scalar>::max();
  std::size_t ind = points.size();

  for (std::size_t i=0; i<points.size(); ++i)
  {
    Scalar dist = points[i]->dot(direction);
    if (dist > max_dist)
    {
      max_dist = dist;
      ind = i;
    }
  }
  std::cout << "mdist = " << max_dist << std::endl;

  assert( ind < points.size() );
  return *points[ind];
}


template<typename Scalar>
std::vector<Point<3,Scalar> *> &
Shape<Scalar>::get_points()
{
  return points;
}

}
