// Geometrical plane object in 3D
#pragma once

#include "Point.hpp"


namespace angem
{

template <int dim, typename Scalar>
class Line
{
 public:
  Line();
  Line(const Point<dim,Scalar> & point,
       const Point<dim,Scalar> & direction);

  Point<dim, Scalar> point;
  Point<dim, Scalar> direction;
};


template <int dim, typename Scalar>
Line<dim,Scalar>::Line()
{}


template <int dim, typename Scalar>
Line<dim,Scalar>::Line(const Point<dim,Scalar> & point,
                       const Point<dim,Scalar> & direction)
    :
    point(point),
    direction(direction)
{}

}  // end namespace
