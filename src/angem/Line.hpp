// Geometrical plane object in 3D
#pragma once

#include "Point.hpp"


namespace angem
{

/* Simple utility class defining a Line object.
 * A line is assigned by a direction vector and and point on the line.
 * I think this should rather be a struct.
 */
template <int dim, typename Scalar>
class Line
{
 public:
  // Default constructor. creates an invalid line. Sometimes useful.
  Line();
  // Creates a line from a point and direction vectors
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
