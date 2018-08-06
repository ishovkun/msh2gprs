#pragma once
#include "Shape.hpp"

namespace angem
{

template<typename Scalar>
class Rectangle: public Shape<Scalar>
{
public:
	Rectangle(Point<3,Scalar> center,
            Scalar          length,
            Scalar          height,
            Scalar          dip_angle,
            Scalar          strike_angle);

 protected:
  std::vector<Point<3,Scalar>> v_points;
  // Plane<Scalar> plane;
};


template<typename Scalar>
Rectangle<Scalar>::Rectangle(Point<3,Scalar> center,
                             Scalar          length,
                             Scalar          height,
                             Scalar          dip_angle,
                             Scalar          strike_angle)
    :
    Shape<Scalar>::Shape()
{
  // convert to radians
  Scalar rdip    = dip_angle * M_PI / 180.;
  Scalar rstrike = strike_angle * M_PI / 180.;

  // define two unit vectors within the square plane
  // heading strike
  Point<3,Scalar> t1 = {cos(rdip), sin(rdip), 0};
  // heading up the fracture
  Point<3,Scalar> t2 = {sin(rstrike), -cos(rdip)*cos(rstrike), sin(rdip)};

  // define rectangle vertices
  // top left
  v_points.emplace_back();
  v_points.back() =
      center - 0.5*length*t1 + 0.5*height*t2;
  // bottom left
  v_points.emplace_back();
  v_points.back() =
      center - 0.5*length*t1 - 0.5*height*t2;
  // bottom right
  v_points.emplace_back();
  v_points.back() =
      center + 0.5*length*t1 - 0.5*height*t2;
  // top right
  v_points.emplace_back();
  v_points.back() =
      center + 0.5*length*t1 + 0.5*height*t2;

  this->set_data(v_points);
  // plane.set_data(v_points[0], v_points[1], v_points[2]);
}

} // end angem
