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
            Scalar          width,
            Scalar          dip_angle,
            Scalar          strike_angle);
};


template<typename Scalar>
Rectangle<Scalar>::Rectangle(Point<3,Scalar> center,
                             Scalar          length,  // along dip
                             Scalar          width,   // in perp direction
                             Scalar          dip_angle,
                             Scalar          strike_angle)
    :
    Polygon<Scalar>::Polygon()
{
  assert(length > 0);
  assert(width > 0);
  // convert to radians
  Scalar rdip    = dip_angle * M_PI / 180.;
  Scalar rstrike = strike_angle * M_PI / 180.;

  // define two unit vectors within the square plane
  // heading strike
  Point<3,Scalar> t1 = {cos(rstrike), sin(rstrike), 0};
  // heading up the fracture (in direction opposit to dipping)
  Point<3,Scalar> t2 = {cos(rdip)*sin(rstrike), -cos(rdip)*cos(rstrike), sin(rdip)};

  std::vector<Point<3,Scalar>> v_points;
  // define rectangle vertices
  // top left
  v_points.emplace_back();
  v_points.back() = center - 0.5*length*t1 + 0.5*width*t2;
  // bottom left
  v_points.emplace_back();
  v_points.back() = center - 0.5*length*t1 - 0.5*width*t2;
  // bottom right
  v_points.emplace_back();
  v_points.back() = center + 0.5*length*t1 - 0.5*width*t2;
  // top right
  v_points.emplace_back();
  v_points.back() = center + 0.5*length*t1 + 0.5*width*t2;

  Polygon<Scalar>::set_data(v_points);

  Basis<3,Scalar> basis({-t2, t1, cross(t1, t2)});
  this->plane.set_basis(basis);
}

} // end angem
