// Generic dimension-independent implementation of a point class

#pragma once

#include <ostream>
#include <cmath>
#include <vector>
#include <cassert>

namespace angem
{

template<int dim=3, typename Scalar=double>
class Point
{
 public:
  // constructors
  Point();
  Point(const Point<dim, Scalar> & p);
  Point(const std::vector<Scalar> & v);
  Point(const Scalar x, const Scalar y);                  // 2D only
  Point(const Scalar x, const Scalar y, const Scalar z);  // 3D only
  // This also works !!!
  // Point<3, double> p = {1, 2, 3}

  // assign
  void operator=(std::vector<Scalar> & v);

  // getters (non-const)
  Scalar & x();
  Scalar & y();
  Scalar & z();
  Scalar & operator[] (int i);

  // getters (const)
  Scalar x() const;
  Scalar y() const;
  Scalar z() const;
  Scalar operator() (int i) const;

  // operations
  // comparison
  bool operator==(const Point<dim, Scalar> & p) const;
  bool operator!=(const Point<dim, Scalar> & p) const;
  // point-wise sum
  void operator+=(const Point<dim, Scalar> & p);
  // point-wise difference
  void operator-=(const Point<dim, Scalar> & p);
  // add to each component
  void operator+=(const Scalar x);
  // subtract from each component
  void operator-=(const Scalar x);
  // component-wise multiplication
  void operator*=(const Scalar x);
  // component-wise division
  void operator/=(const Scalar x);
  // dot product
  Scalar dot(const Point<dim, Scalar> & p) const;
  // cross product -- only 3D
  void cross(const Point<3, Scalar> & p,
             Point<3, Scalar>       & result) const;
  Point<3, Scalar> cross(const Point<3, Scalar> & p) const;

  // Euclidian distance
  Scalar distance(const Point<dim, Scalar> & p) const;
  // Euclidian norm
  Scalar norm() const;
  // divide by norm
  void normalize();

  // also this class defines the following external functions
  // point-wise sum
  // Point<dim, Scalar> operator+(const Point<dim, Scalar> & p1,
  //                              const Point<dim, Scalar> & p2);

  // point-wise subtraction
  // Point<dim, Scalar> operator-(const Point<dim, Scalar> & p1,
  //                              const Point<dim, Scalar> & p2);

  // dot product
  // Scalar operator*(const Point<dim, Scalar> & p1,
  //                  const Point<dim, Scalar> & p2);

  // Euclidian distance
  // Scalar distance(const Point<dim, Scalar> & p1,
  //                 const Point<dim, Scalar> & p2);

  // Euclidian norm
  // Scalar norm(const Point<dim, Scalar> & p1);

  // printout
  // template<int dim, typename Scalar>
  // std::ostream &operator<<(std::ostream            & os,
  //                          const Point<dim,Scalar> & p)


 protected:
  Scalar coords[dim];
};


// CONSTRUCTORS
template<int dim, typename Scalar>
Point<dim,Scalar>::Point()
{
  assert(dim > 0 and dim <= 3);
  for (int i=0; i<dim; ++i)
    coords[i] = 0;
}


template<int dim, typename Scalar>
Point<dim,Scalar>::Point(const Point<dim, Scalar> & p)
{
  assert(dim > 0 and dim <= 3);
  for (int i=0; i<dim; ++i)
    coords[i] = p.coords[i];
}


template<int dim, typename Scalar>
Point<dim,Scalar>::Point(const std::vector<Scalar> & v)
{
  assert(dim > 0 and dim <= 3);
  assert(v.size() == dim);
  for (int i=0; i<dim; ++i)
    coords[i] = v[i];
}


// Assign
template<int dim, typename Scalar>
void Point<dim,Scalar>::operator=(std::vector<Scalar> & v)
{
  assert(v.size() == dim);
  for (int i=0; i<dim; ++i)
    coords[i] = v[i];
}


// Partial specialization -- 2D
template<int dim, typename Scalar>
Point<dim,Scalar>::Point(const Scalar x, const Scalar y)
{
  assert(dim == 2);
  this->coords[0] = x;
  this->coords[1] = y;
}


template<int dim, typename Scalar>
Point<dim,Scalar>::Point(const Scalar x, const Scalar y, const Scalar z)
{
  assert(dim == 3);
  coords[0] = x;
  coords[1] = y;
  coords[2] = z;
}


// GETTERS (CONST)
template<int dim, typename Scalar>
Scalar Point<dim,Scalar>::x() const
{
  return coords[0];
}


template<int dim, typename Scalar>
Scalar Point<dim,Scalar>::y() const
{
  assert(dim > 1);
  return coords[1];
}


template<int dim, typename Scalar>
Scalar Point<dim,Scalar>::z() const
{
  assert(dim == 3);
  return coords[2];
}


template<int dim, typename Scalar>
Scalar Point<dim,Scalar>::operator() (int i) const
{
  assert(i < dim);
  return coords[i];
}


// GETTERS (non-CONST)
template<int dim, typename Scalar>
Scalar & Point<dim,Scalar>::x()
{
  return coords[0];
}


template<int dim, typename Scalar>
Scalar & Point<dim,Scalar>::y()
{
  assert(dim > 1);
  return coords[1];
}


template<int dim, typename Scalar>
Scalar & Point<dim,Scalar>::z()
{
  assert(dim == 3);
  return coords[2];
}


template<int dim, typename Scalar>
Scalar & Point<dim,Scalar>::operator[] (int i)
{
  assert(i < dim);
  return coords[i];
}


// OPERATORS
template<int dim, typename Scalar>
bool Point<dim,Scalar>::operator==(const Point<dim, Scalar> & p) const
{
  for (int i=0; i<dim; ++i)
  {
    if (coords[i] == p(i))
      continue;
    else
      return false;
  }
  return true;
}


template<int dim, typename Scalar>
bool Point<dim,Scalar>::operator!=(const Point<dim, Scalar> & p) const
{
  for (int i=0; i<dim; ++i)
  {
    if (coords[i] == p(i))
      continue;
    else
      return true;
  }
  return false;
}


template<int dim, typename Scalar>
void Point<dim,Scalar>::operator+=(const Point<dim, Scalar> & p)
{
  for (int i=0; i<dim; ++i)
    coords[i] += p(i);
}


template<int dim, typename Scalar>
void Point<dim,Scalar>::operator-=(const Point<dim, Scalar> & p)
{
  for (int i=0; i<dim; ++i)
    coords[i] -= p(i);
}


template<int dim, typename Scalar>
void Point<dim,Scalar>::operator+=(const Scalar x)
{
  for (int i=0; i<dim; ++i)
    coords[i] += x;
}


template<int dim, typename Scalar>
void Point<dim,Scalar>::operator-=(const Scalar x)
{
  for (int i=0; i<dim; ++i)
    coords[i] -= x;
}


template<int dim, typename Scalar>
void Point<dim,Scalar>::operator*=(const Scalar x)
{
  for (int i=0; i<dim; ++i)
    coords[i] *= x;
}


template<int dim, typename Scalar>
void Point<dim,Scalar>::operator/=(const Scalar x)
{
  for (int i=0; i<dim; ++i)
    coords[i] /= x;
}


template<int dim, typename Scalar>
Scalar Point<dim,Scalar>::dot(const Point<dim, Scalar> & p) const
{
  Scalar result = static_cast<Scalar>(0);
  for (int i=0; i<dim; ++i)
    result += coords[i] * p(i);
  return result;
}


template<int dim, typename Scalar>
void Point<dim,Scalar>::cross(const Point<3, Scalar> & p,
                              Point<3, Scalar>       & result) const
{
  result[0] = coords[1]*p(2) - coords[2]*p(1);
  result[1] = coords[2]*p(0) - coords[0]*p(2);
  result[2] = coords[0]*p(1) - coords[1]*p(0);
}


template<int dim, typename Scalar>
Point<3,Scalar> Point<dim,Scalar>::cross(const Point<3, Scalar> & p) const
{
  Point<3,Scalar> result;
  cross(p, result);
  return result;
}


template<int dim, typename Scalar>
Scalar Point<dim,Scalar>::distance(const Point<dim, Scalar> & p) const
{
  Scalar result = static_cast<Scalar>(0);
  for (int i=0; i<dim; ++i)
  {
    const Scalar d = coords[i] - p(i);
    result += d * d;
  }
  return std::sqrt(result);
}


template<int dim, typename Scalar>
Scalar Point<dim,Scalar>::norm() const
{
  Scalar result = static_cast<Scalar>(0);
  for (int i=0; i<dim; ++i)
  {
    const Scalar d = coords[i];
    result += d * d;
  }
  return std::sqrt(result);
}


template<int dim, typename Scalar>
void Point<dim,Scalar>::normalize()
{
  Scalar nor = norm();
  assert(nor != static_cast<Scalar>(0));
  for (int i=0; i<dim; ++i)
    coords[i] /= nor;
}

// EXTERNAL OPERATORS
template<int dim, typename Scalar>
Point<dim, Scalar> operator+(const Point<dim, Scalar> & p1,
                             const Point<dim, Scalar> & p2)
{
  Point<dim, Scalar> result;
  for (int i=0; i<dim; ++i)
    result[i] = p1(i) + p2(i);
  return result;
}


template<int dim, typename Scalar>
Point<dim, Scalar> operator-(const Point<dim, Scalar> & p1,
                             const Point<dim, Scalar> & p2)
{
  Point<dim, Scalar> result;
  for (int i=0; i<dim; ++i)
    result[i] = p1(i) - p2(i);
  return result;
}


template<int dim, typename Scalar>
Point<dim, Scalar> operator-(const Point<dim, Scalar> & p)
{
  Point<dim, Scalar> result;
  for (int i=0; i<dim; ++i)
    result[i] = -p(i);
  return result;
}


template<int dim, typename Scalar>
Scalar operator*(const Point<dim, Scalar> & p1,
                 const Point<dim, Scalar> & p2)
{
  return p1.dot(p2);
}


template<int dim, typename Scalar>
Point<dim, Scalar> operator*(const Point<dim, Scalar> & p,
                             const Scalar             & x)
{
  Point<dim, Scalar> result;
  for (int i=0; i<dim; ++i)
    result[i] = p(i) * x;
  return result;
}


template<int dim, typename Scalar>
Point<dim, Scalar> operator*(const Scalar             & x,
                             const Point<dim, Scalar> & p)
{
  return p*x;
}


template<int dim, typename Scalar>
Point<dim, Scalar> operator/(const Point<dim, Scalar> & p,
                             const Scalar             & x)
{
  Point<dim, Scalar> result;
  for (int i=0; i<dim; ++i)
    result[i] = p(i) / x;
  return result;
}


template<int dim, typename Scalar>
Scalar distance(const Point<dim, Scalar> & p1,
                const Point<dim, Scalar> & p2)
{
  return p1.distance(p2);
}


template<int dim, typename Scalar>
Scalar norm(const Point<dim, Scalar> & p1)
{
  return p1.norm();
}


template<typename Scalar>
void cross_product(const Point<3, Scalar> & p1,
                   const Point<3, Scalar> & p2,
                   Point<3, Scalar>       & result)
{
  p1.cross(p2, result);
}


template<typename Scalar>
Point<3,Scalar> cross_product(const Point<3, Scalar> & p1,
                              const Point<3, Scalar> & p2)
{
  Point<3,Scalar> result;
  p1.cross(p2, result);
  return result;
}

// PRINTING
template<int dim, typename Scalar>
std::ostream &operator<<(std::ostream            & os,
                         const Point<dim,Scalar> & p)
{
  for (int i=0; i<dim; ++i)
    os << p(i) << "\t";
  return os;
}

}  // end namespace
