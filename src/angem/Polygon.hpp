#pragma once

#include <Shape.hpp>
#include <Plane.hpp>
#include <utils.hpp>

namespace angem
{

template<typename Scalar>
class Polygon: public Shape<Scalar>
{
 public:
	Polygon();
  Polygon(std::vector<Point<3,Scalar> *> & points_list);
  Point<3,Scalar> center() const {return plane.point;};

  // shift all points in direction p
  virtual void set_data(std::vector<Point<3,Scalar>> & point_list) override;
  virtual void move(const Point<3,Scalar> & p) override;

  // Attributes
  Plane<Scalar> plane;
};


template<typename Scalar>
Polygon<Scalar>::Polygon()
    :
    Shape<Scalar>::Shape()
{}


template<typename Scalar>
Polygon<Scalar>::Polygon(std::vector<Point<3,Scalar> *> & point_list)
    :
    Shape<Scalar>::Shape(point_list)
{}


template<typename Scalar>
void
Polygon<Scalar>::set_data(std::vector<Point<3,Scalar>> & point_list)
{
  Shape<Scalar>::set_data(point_list);

  assert(point_list.size() > 3);
  Point<3,Scalar> center = compute_center_mass(point_list);
  plane.set_data(compute_center_mass(point_list),
                 point_list[1], point_list[2]);
}


template<typename Scalar>
void
Polygon<Scalar>::move(const Point<3,Scalar> & p)
{
  Shape<Scalar>::move(p);
  plane.move(p);
}
}
