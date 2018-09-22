// Geometrical plane object in 3D
#pragma once

#include "Point.hpp"
#include "Basis.hpp"
#include "Line.hpp"

#include <math.h>    // M_PI
#include <algorithm> // clamp

namespace angem
{

template <typename Scalar>
class Plane
{
 public:
  Plane();
  // create plane from a point on the plane and a normal vector
  Plane(const Point<3,Scalar> & point,
        const Point<3,Scalar> & normal);
  // strike and dip in degrees
  Plane(const Point<3,Scalar> & point,
        const Scalar          & dip_angle,      // -90 <= dip <= 90
        const Scalar          & strike_angle);
  // create plane from 3 points
  Plane(const Point<3,Scalar> & p1,
        const Point<3,Scalar> & p2,
        const Point<3,Scalar> & p3);

  // setter
  void set_data(const Point<3,Scalar> & p1,
                const Point<3,Scalar> & p2,
                const Point<3,Scalar> & p3);

  // shift support in direction p
  void move(const Point<3,Scalar> & p);

  const Point<3,Scalar> & normal() const {return basis(2);}
  Basis<3,Scalar> & get_basis() {return basis;}

  // compute strike angle (degrees) from normal
  Scalar strike_angle() const;
  // compute dip angle (degrees) from normal
  Scalar dip_angle() const;

  // get point projection on the plane (in global coordinates)
  Point<3,Scalar> project(const Point<3,Scalar> & p) const;
  Point<3,Scalar> local_coordinates(const Point<3,Scalar> & p) const;
  // project vector (no account for plane location)
  Point<3,Scalar> project_vector(const Point<3,Scalar> & p) const;

  // signed distance from point to plane (> 0 if point is above plane)
  Scalar distance(const Point<3,Scalar> & p) const;

  // true if point is above the plane
  bool above(const Point<3,Scalar> & p) const;

  // ATTRIBUTES
  Point<3,Scalar> point;  // point on the plane
  // algebraic coefficient in ax + by + cz = d
  Scalar d;

 protected:
  void compute_algebraic_coeff();
  // return two orthogonal vectors within the plane
  void compute_basis(const Point<3,Scalar> & normal);

  Basis<3,Scalar> basis;
};


template <typename Scalar>
Plane<Scalar>::Plane()
{}


template <typename Scalar>
Plane<Scalar>::Plane(const Point<3,Scalar> & point,
                     const Point<3,Scalar> & normal)
    :
    point(point)
{
  assert(fabs(normal.norm() - 1) < 1e-12);
  compute_basis(normal);
  compute_algebraic_coeff();
}


template <typename Scalar>
Plane<Scalar>::Plane(const Point<3,Scalar> & point,
                     const Scalar          & dip_angle,
                     const Scalar          & strike_angle)
    :
    point(point)
{
  assert(dip_angle >= - 90 and dip_angle <= 90);

  Scalar rdip    = dip_angle * M_PI / 180.;
  Scalar rstrike = strike_angle * M_PI / 180.;

  Point<3,Scalar> normal;
  // dip = polar angle between normal and z
  normal[0] = sin(rdip) * cos(rstrike + M_PI/2.);
  normal[1] = sin(rdip) * sin(rstrike + M_PI/2.);
  normal[2] = cos(rdip);

  compute_basis(normal);
  compute_algebraic_coeff();
}


template <typename Scalar>
Plane<Scalar>::Plane(const Point<3,Scalar> & p1,
                     const Point<3,Scalar> & p2,
                     const Point<3,Scalar> & p3)
{
  set_data(p1, p2, p3);
}


template <typename Scalar>
void Plane<Scalar>::set_data(const Point<3,Scalar> & p1,
                             const Point<3,Scalar> & p2,
                             const Point<3,Scalar> & p3)
{
  assert(p1 != p2 and p2 != p3 and p1 != p3);

  point = p1;
  // define two tangent vectors
  const Point<3> t1 = p1 - p2;
  const Point<3> t2 = p1 - p3;
  Point<3,Scalar> normal = cross(t1, t2);
  normal.normalize();

  compute_basis(normal);
  compute_algebraic_coeff();
}


template <typename Scalar>
void Plane<Scalar>::compute_algebraic_coeff()
{
  const auto & normal = basis[2];
  d = normal(0)*point(0) +
      normal(1)*point(1) +
      normal(2)*point(2);
}


template <typename Scalar>
Scalar Plane<Scalar>::distance(const Point<3,Scalar> & p) const
{
  /* dot product of point by perpendicular vector:
   * if < 0: point is below the plane
   */
  /* signed distance from point p (vertex) to plane
   * with normal n and containing point x0:
   * d = (p - x0) · n
   * if d<0, point below the plane
   */
  return (p - point).dot(normal());
}


template <typename Scalar>
bool Plane<Scalar>::above(const Point<3,Scalar> & p) const
{
  if (distance(p) > 0)
    return true;
  else
    return false;
}


// Determine if a set of points lies within a plane
template <typename Scalar>
bool align_on_plane(const std::vector<Point<3, Scalar>> & points)
{
  assert(points.size() > 2);

  if (points.size() == 3)
    return true;

  Plane<Scalar> plane(points[0], points[1], points[2]);
  for (std::size_t i=3; i<points.size(); ++i)
    if ( fabs(plane.distance(points[i])) > 1e-16 )
      return false;

  return true;
}


// return two orthogonal vectors within the plane
template <typename Scalar>
void Plane<Scalar>::compute_basis(const Point<3,Scalar> & normal)
{
  /*
   * Algorithm:
   * chose a random vector rv
   * project it on the plane - first tangent basis vector e1
   * complete the basis with e2 = n x e1
   */
  // const Point<3,Scalar> rv = {normal.x() + 1, normal.y(), normal.z()};
  Point<3,Scalar> rv = normal;
  rv[0] += 1;
  Point<3,Scalar> e1 = project_vector(rv);

  // check
  // assert(e1.norm() > 1e-16);
  if (e1.norm() < 1e-16)
  {
    rv = normal;
    rv[1] += 1;
    e1 = project_vector(rv);
  }
  e1.normalize();

  Point<3,Scalar> e2 = normal.cross(e1);
  e2.normalize();

  basis[0] = e1;
  basis[1] = e2;
  basis[2] = normal;
}


// project the direction of the vector on the plane (no translation)
template <typename Scalar>
inline
Point<3, Scalar>
Plane<Scalar>::project_vector(const Point<3,Scalar> & p) const
{
  Scalar p_n = p.dot(basis(2));
  return p - p_n * basis(2);
}


// get point projection on the plane (in global coordinates considering plane translation)
template <typename Scalar>
inline
Point<3, Scalar>
Plane<Scalar>::project(const Point<3,Scalar> & p) const
{
  // 1: translate p' = p - s  (s - plane support point)
  // 2. project on normal p'n = p' · n
  // 3. Translate back j = s + j'
  Point<3,Scalar> p_prime = p - point;
  Point<3,Scalar> j_prime = project_vector(p_prime);
  return point + j_prime;
}


// get point projection on the plane (in both global coordinates and local basis)
template <typename Scalar>
Point<3,Scalar>
Plane<Scalar>::local_coordinates(const Point<3,Scalar> & p) const
{
  // translate
  Point<3,Scalar> p_prime = p - point;
  // project on basis vectors
  return basis.transform(p_prime);
}


template <typename Scalar>
void
Plane<Scalar>::move(const Point<3,Scalar> & p)
{
  point += p;
  compute_algebraic_coeff();
}


template <typename Scalar>
Scalar
Plane<Scalar>::dip_angle() const
{
  Scalar rdip = static_cast<Scalar>(acos(basis(2)[2]));
  double dip = 180. * rdip / M_PI;
  if (dip > 90.0)
    dip = 180. - dip;
  return dip;
}


template <typename Scalar>
Scalar
Plane<Scalar>::strike_angle() const
{
  Scalar rdip = static_cast<Scalar>(acos(basis(2)[2]));

  // avoid taking acos(+- 1) -- causes errors due to roundoff
  const double v1 = std::clamp( basis(2)[0] / sin(rdip), -1.0, 1.0);
  const double v2 = std::clamp( basis(2)[1] / sin(rdip), -1.0, 1.0);

  Scalar rstrike_from_cos = acos(v1) - M_PI / 2.;
  Scalar rstrike_from_sin = asin(v2) - M_PI / 2.;
  std::cout << "rstrike_from_sin = "<< rstrike_from_sin << std::endl;
  std::cout << "rstrike_from_cos = "<< rstrike_from_cos << std::endl;


  Scalar strike;
  if (rstrike_from_sin >= 0 and rstrike_from_cos >= 0)
  {
    strike = 180. * rstrike_from_cos / M_PI;
  }
  else if (rstrike_from_sin < 0 and rstrike_from_cos > 0)
  {
    strike = - 180. * rstrike_from_cos / M_PI;
  }
  else if (rstrike_from_sin > 0 and rstrike_from_cos < 0)
  {
    strike = 180. * rstrike_from_cos / M_PI;
    strike = fabs(strike);
  }
  else // if (rstrike_from_sin < 0 and rstrike_from_cos < 0)
  {
    strike = 180. * rstrike_from_cos / M_PI;
    strike = fabs(strike);
  }

  return strike;
}
}  // end namespace
