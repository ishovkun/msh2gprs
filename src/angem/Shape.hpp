#pragma once

#include <Point.hpp>
#include <limits>

namespace angem
{

template<typename Scalar>
class Shape
  {
   public:
    // constructors
    Shape();
    Shape(std::vector<Point<3,Scalar>> & point_list);
    Shape(std::vector<Point<3,Scalar> *> & points_list);
    // setters
    virtual void set_data(std::vector<Point<3,Scalar>> & point_list);
    virtual void set_data(std::vector<Point<3,Scalar> *> & point_list);

    // getter
    std::vector<Point<3,Scalar> *> & get_points();
    // check if empty
    bool empty() const;
    // support function for gjk collision algorithm
    virtual Point<3,Scalar> support(Point<3, Scalar> & direction) const;
    // shift all points in direction p
    virtual void move(const Point<3,Scalar> & p);

   protected:
    std::vector<Point<3,Scalar> *> points;
};


template<typename Scalar>
Shape<Scalar>::Shape()
{}


template<typename Scalar>
Shape<Scalar>::Shape(std::vector<Point<3,Scalar> *> & point_list)
    :
    points(point_list)
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
void
Shape<Scalar>::set_data(std::vector<Point<3,Scalar> *> & point_list)
{
  points = point_list;
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

  assert( ind < points.size() );
  return *points[ind];
}


template<typename Scalar>
std::vector<Point<3,Scalar> *> &
Shape<Scalar>::get_points()
{
  return points;
}


template<typename Scalar>
void
Shape<Scalar>::move(const Point<3,Scalar> & p)
{
  for (std::size_t i=0; i<points.size(); ++i)
    *points[i] += p;
}

}
