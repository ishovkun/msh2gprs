#pragma once

#include <Point.hpp>
#include <limits>
#include <iostream>

namespace angem
{

template<typename Scalar>
class Shape
  {
   public:
    // constructors
    Shape();
    Shape(const std::vector<const Point<3,Scalar>> & point_list);
    Shape(const std::vector<Point<3,Scalar>> & all_mesh_vertices,
          const std::vector<std::size_t>     & indices);
    // setters
    virtual void set_data(const std::vector<Point<3,Scalar>> & point_list);

    // getter
    std::vector<Point<3,Scalar>> & get_points();
    const std::vector<Point<3,Scalar>> & get_points() const;
    // check if empty
    bool empty() const;
    // support function for gjk collision algorithm
    virtual Point<3,Scalar> support(const Point<3, Scalar> & direction) const;
    // center mass
    virtual Point<3,Scalar> center() const;
    // shift all points in direction p
    virtual void move(const Point<3,Scalar> & p);

   protected:
    std::vector<Point<3,Scalar>> points;
};


template<typename Scalar>
Shape<Scalar>::Shape()
{}


template<typename Scalar>
Shape<Scalar>::Shape(const std::vector<const Point<3,Scalar>> & point_list)
    :
    points(point_list)
{}


template<typename Scalar>
void
Shape<Scalar>::set_data(const std::vector<Point<3,Scalar>> & point_list)
{
  points = point_list;
}


template<typename Scalar>
Shape<Scalar>::Shape(const std::vector<Point<3,Scalar>> & all_mesh_vertices,
                     const std::vector<std::size_t>     & indices)
{
  for (const std::size_t & i : indices)
    points.push_back(all_mesh_vertices[i]);
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
Shape<Scalar>::support(const Point<3,Scalar> & direction) const
{
  assert(!empty());

  // get the farthest point in direction p
  Scalar max_dist = - std::numeric_limits<Scalar>::max();
  std::size_t ind = points.size();

  for (std::size_t i=0; i<points.size(); ++i)
  {
    Scalar dist = points[i].dot(direction);
    if (dist > max_dist)
    {
      max_dist = dist;
      ind = i;
    }
  }

  assert( ind < points.size() );
  return points[ind];
}


template<typename Scalar>
std::vector<Point<3,Scalar>> &
Shape<Scalar>::get_points()
{
  return points;
}


template<typename Scalar>
const std::vector<Point<3,Scalar>> &
Shape<Scalar>::get_points() const
{
  return points;
}


template<typename Scalar>
void
Shape<Scalar>::move(const Point<3,Scalar> & p)
{
  for (std::size_t i=0; i<points.size(); ++i)
    points[i] += p;
}


template<typename Scalar>
Point<3,Scalar>
Shape<Scalar>::center() const
{
  return compute_center_mass(points);
}


template<typename Scalar>
std::ostream &operator<<(std::ostream        & os,
                         const Shape<Scalar> & p)
{
  const std::vector<Point<3,Scalar>> & points = p.get_points();

  for (std::size_t i=0; i<points.size(); ++i)
  {
    os << points[i];
    if (i < points.size() - 1)
      os << std::endl;
  }
  return os;
}

}
