#pragma once
#include <Polygon.hpp>

namespace angem
{

template<typename Scalar>
class Rectangle: public Polygon<Scalar>
{
public:
	Rectangle(Point<3,Scalar> center,
            Scalar          length,
            Scalar          height,
            Scalar          dip_angle,
            Scalar          strike_angle);

 protected:
  std::vector<Point<3,Scalar>> v_points;
};


template<typename Scalar>
Rectangle<Scalar>::Rectangle(Point<3,Scalar> center,
                             Scalar          length,
                             Scalar          height,
                             Scalar          dip_angle,
                             Scalar          strike_angle)
    :
    Polygon<Scalar>::Polygon()
{
  assert(length > 0);
  assert(height > 0);
  // convert to radians
  Scalar rdip    = dip_angle * M_PI / 180.;
  Scalar rstrike = strike_angle * M_PI / 180.;

  // define two unit vectors within the square plane
  // heading strike
  Point<3,Scalar> t1 = {cos(rstrike), sin(rstrike), 0};
  // heading up the fracture
  Point<3,Scalar> t2 = {cos(rdip)*sin(rstrike), -cos(rdip)*cos(rstrike), sin(rdip)};

  v_points.clear();
  // define rectangle vertices
  // top left
  v_points.emplace_back();
  v_points.back() = center - 0.5*length*t1 + 0.5*height*t2;
  // bottom left
  v_points.emplace_back();
  v_points.back() = center - 0.5*length*t1 - 0.5*height*t2;
  // bottom right
  v_points.emplace_back();
  v_points.back() = center + 0.5*length*t1 - 0.5*height*t2;
  // top right
  v_points.emplace_back();
  v_points.back() = center + 0.5*length*t1 + 0.5*height*t2;

  Polygon<Scalar>::set_data(v_points);
}

} // end angem
