// Generic dimension-independent implementation of a point class

#pragma once

#include <iostream>  // debug
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
  // CONSTRUCTORS
  // Default constructor. Creates a point with all components equal to zero
  Point();
  // Copy constructor
  Point(const Point<dim,Scalar> & p);
  // Convenience constructor.
  // Creates a point from std::vector. Asserts that vector.size() == dim
  Point(const std::vector<Scalar> & v);
  // Constructor for a 2D Point (dim=2) only.
  // Creates a point with the coordinates equal to x and y, respectively.
  Point(const Scalar x, const Scalar y);                  // 2D only
  // Constructor for a 3D Point (dim=3) only.
  // Creates a point with the coordinates equal to x, y, and z, respectively.
  Point(const Scalar x, const Scalar y, const Scalar z);  // 3D only

  // assignment operator
  void operator=(std::vector<Scalar> & v);
  // set all component to zero
  void clear();

  // getters (non-const)
  // Get a non-const reference to x component
  Scalar & x();
  // Get a non-const reference to y component (dim=2 || dim=3 only)
  Scalar & y();
  // Get a non-const reference to z component (dim=3 only)
  Scalar & z();
  // Get a non-const reference to  i-th coordinate. checks that i < dim
  Scalar & operator[] (int i);

  // getters (const)
  // Get a const reference to x component
  Scalar x() const;
  // Get a const reference to y component (dim=2 || dim=3 only)
  Scalar y() const;
  // Get a const reference to z component (dim=3 only)
  Scalar z() const;
  // Get a const reference to  i-th coordinate. checks that i < dim
  Scalar operator() (int i) const;
  // Get a const reference to  i-th coordinate. checks that i < dim
  const Scalar & operator[] (int i) const;

  // operations
  // comparison
  // returns true if all components are the same
  bool operator==(const Point<dim, Scalar> & p) const;
  // returns true if nay components are different
  bool operator!=(const Point<dim, Scalar> & p) const;
  // compares the norms, can be used for sorting
  bool operator< (const Point<dim, Scalar> & other) const;
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

  // if norm of the cross_product is small
  bool parallel(const Point<dim, Scalar> & other,
                const double               tol = 1e-5) const;

  // Euclidian distance
  Scalar distance(const Point<dim, Scalar> & p) const;
  // Euclidian norm
  Scalar norm() const;
  // divide by norm
  void normalize();

  // also this class defines the following external functions
  // point-wise sum
  template <int d, typename S>
  friend Point<d, S> operator+(const Point<d,S> & p1,
                               const Point<d,S> & p2);
  // point-wise subtraction
  template <int d, typename S>
  friend Point<d, S> operator-(const Point<d,S> & p1,
                               const Point<d,S> & p2);
  // multiplication by number
  template <int d, typename S1, typename S2>
  friend Point<d, S1> operator*(const Point<d,S1> & p1,
                                const S2          & x);
  // division by number
  template <int d, typename S1, typename S2>
  friend Point<d, S1> operator/(const Point<d,S1> & p1,
                                const S2          & x);
  // dot product
  template <int d, typename S>
  friend S operator*(const Point<d,S> & p1,
                     const Point<d,S> & p2);

  // Euclidian distance
  template <int d, typename S>
  friend S distance(const Point<d,S> & p1,
                    const Point<d,S> & p2);

  // Euclidian norm
  template <int d, typename S>
  friend S norm(const Point<d,S> & p1);

  // printout
  template <int d, typename S>
  friend std::ostream &operator<<(std::ostream     & os,
                                  const Point<d,S> & p);


 protected:
  Scalar coords[dim];  // array of coordinate components
};


// CONSTRUCTORS
template<int dim, typename Scalar>
Point<dim,Scalar>::Point()
{
  static_assert(dim > 0 and dim <= 3,
                "Only 1,2, and 3 dimensions are supported");
  for (int i=0; i<dim; ++i)
    coords[i] = 0;
}


template<int dim, typename Scalar>
Point<dim,Scalar>::Point(const Point<dim, Scalar> & p)
{
  static_assert(dim > 0 and dim <= 3,
                "Only 1,2, and 3 dimensions are supported");
  for (int i=0; i<dim; ++i)
    coords[i] = p.coords[i];
}


template<int dim, typename Scalar>
Point<dim,Scalar>::Point(const std::vector<Scalar> & v)
{
  static_assert(dim > 0 and dim <= 3,
                "Only 1,2, and 3 dimensions are supported");
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


template<int dim, typename Scalar>
void Point<dim,Scalar>::clear()
{
  for (int i=0; i<dim; ++i)
    coords[i] = 0;
}


// Partial specialization -- 2D
template<int dim, typename Scalar>
Point<dim,Scalar>::Point(const Scalar x, const Scalar y)
{
  static_assert(dim == 2, "Only 2d objects can be initialized this way");
  this->coords[0] = x;
  this->coords[1] = y;
}


template<int dim, typename Scalar>
Point<dim,Scalar>::Point(const Scalar x, const Scalar y, const Scalar z)
{
  static_assert(dim == 3,
                "Only 3d objects can be initialized this way");
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
  static_assert(dim > 1,
                "1d objects have only one coordinate");
  return coords[1];
}


template<int dim, typename Scalar>
Scalar Point<dim,Scalar>::z() const
{
  static_assert(dim == 3,
                "Only 3D objects have z coordinate");
  return coords[2];
}


template<int dim, typename Scalar>
Scalar Point<dim,Scalar>::operator() (int i) const
{
  assert(i < dim);
  return coords[i];
}


template<int dim, typename Scalar>
const Scalar & Point<dim,Scalar>::operator[] (int i) const
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
  static_assert(dim == 3,
                "Only 3D objects have z coordinate");
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
bool Point<dim,Scalar>::operator< (const Point<dim, Scalar> & other) const
{
  if (norm() < other.norm())
    return true;
  else
  {
    return false;
  }
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
bool Point<dim,Scalar>::parallel(const Point<dim, Scalar> & other,
                                 const double               tol) const
{
  // auto result = (cross(other)).norm();
  // std::cout << "res = " << result << std::endl;
  if ( (cross(other)).norm() < tol )
  {
    return true;
  }
  else
    return false;
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
  assert (std::isfinite(static_cast<Scalar>(1) / nor));
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


template<int dim, typename Scalar1, typename Scalar2>
Point<dim, Scalar1> operator*(const Point<dim, Scalar1> & p1,
                              const Scalar2             & x)
{
  Point<dim, Scalar1> result;
  for (int i=0; i<dim; ++i)
    result[i] = p1.coords[i] * static_cast<Scalar1>(x);
  return result;
}


template<int dim, typename Scalar1, typename Scalar2>
Point<dim, Scalar1> operator/(const Point<dim, Scalar1> & p1,
                              const Scalar2             & x)
{
  Point<dim, Scalar1> result;
  for (int i=0; i<dim; ++i)
    result[i] = p1.coords[i] / static_cast<Scalar1>(x);
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


template<int dim, typename Scalar1, typename Scalar2>
Point<dim, Scalar1> operator*(const Scalar2             & x,
                              const Point<dim, Scalar1> & p)
{
  return p * x;
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
Scalar dot(const Point<3, Scalar> & p1,
           const Point<3, Scalar> & p2)
{
  return p1.dot(p2);
}


template<typename Scalar>
Point<3,Scalar> cross(const Point<3, Scalar> & p1,
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
  {
    os << p(i);
    if (i < dim - 1)
      os << "\t";
  }
  return os;
}


template<int dim, typename Scalar>
std::ostream &operator<<(std::ostream                         & os,
                         const std::vector<Point<dim,Scalar>> & points)
{
  for (const auto & p : points)
    os << p << std::endl;
  return os;
}


template<int dim, typename Scalar>
inline
std::size_t
find(const Point<dim,Scalar>              & p,
     const std::vector<Point<dim,Scalar>> & points,
     const double                           tol = 1e-6)
{
  std::size_t counter = 0;
  for (const auto & point : points)
  {
    if (p.distance(point) < tol)
      return counter;
    counter++;
  }
  return counter;
}


template<int dim, typename Scalar>
inline
std::size_t insert(const Point<dim,Scalar>        & point,
                   std::vector<Point<dim,Scalar>> & points,
                   const Scalar                     tol = 1e-6)
{
  const std::size_t ind = find(point, points, tol);
  if (ind == points.size())
    points.push_back(point);
  return ind;
}



}  // end namespace

namespace std
{
template <>
struct hash<angem::Point<3,double>>
{
  size_t operator()(angem::Point<3,double> const & p) const noexcept
  {
    // return (
    //     (51 + std::hash<int>()(p.x())) * 51 +
    //     (23 + std::hash<int>()(p.y())) * 23 +
    //     (std::hash<int>()(p.z()))
    //         );
    return (
        (51 + std::hash<double>()(p.x())) * 51 +
        (23 + std::hash<double>()(p.y())) * 23 +
        (std::hash<double>()(p.z()))
        );
  }
};

}
