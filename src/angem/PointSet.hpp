#pragma once
#include <Point.hpp>

namespace angem
{

/* Set of points
 * Insertion and lookup in linear time since
 * points are compared with a tolerance.
 * I guess it would be cool to implement fast hashing for
 * points with tolerance
 */

template<int dim, typename Scalar>
class PointSet
{
 public:
  PointSet(const double tol = 1e-6);
  // search for point in set and append if not found.
  // returns point index in set
  std::size_t insert(const Point<dim,Scalar> &p);
  // wrapper around find(const Point &p, std::vector<Point> &vPoints)
  std::size_t find(const Point<dim,Scalar> &p);

  // getters
  Point<dim,Scalar> & operator[](std::size_t i);
  const Point<dim,Scalar> & operator[](std::size_t i) const;

  // STD-like operators
  // iterators:
  typename std::vector<Point<dim,Scalar>>::iterator begin();
  typename std::vector<Point<dim,Scalar>>::iterator end();
  typename std::vector<Point<dim,Scalar>>::const_iterator begin() const;
  typename std::vector<Point<dim,Scalar>>::const_iterator end() const;

  // size methods. returns size of vector points
  std::size_t size() const;

  // variables
  std::vector<Point<dim,Scalar>> points;
 private:
  // search tolerance
  double tol;
};


template<int dim, typename Scalar>
PointSet<dim,Scalar>::PointSet(const double tol)
    :
    tol(tol)
{}



template<int dim, typename Scalar>
std::size_t PointSet<dim,Scalar>::find(const Point<dim,Scalar> &p)
{
  return angem::find(p, points, tol);
}


template<int dim, typename Scalar>
std::size_t PointSet<dim,Scalar>::insert(const Point<dim,Scalar> &p)
{
  // return insert(p, points, tol);
  return angem::insert<dim,Scalar>(p, points, tol);
}


template<int dim, typename Scalar>
Point<dim,Scalar> & PointSet<dim,Scalar>::operator[](std::size_t i)
{
  return points[i];
}


template<int dim, typename Scalar>
const Point<dim,Scalar> & PointSet<dim,Scalar>::operator[](std::size_t i) const
{
  return points[i];
}


template<int dim, typename Scalar>
std::size_t PointSet<dim,Scalar>::size() const
{
  return points.size();
}


template<int dim, typename Scalar>
typename std::vector<Point<dim,Scalar>>::iterator PointSet<dim,Scalar>::begin()
{
  return points.begin();
}


template<int dim, typename Scalar>
typename std::vector<Point<dim,Scalar>>::iterator PointSet<dim,Scalar>::end()
{
  return points.end();
}


template<int dim, typename Scalar>
typename std::vector<Point<dim,Scalar>>::const_iterator PointSet<dim,Scalar>::begin() const
{
  return points.begin();
}


template<int dim, typename Scalar>
typename std::vector<Point<dim,Scalar>>::const_iterator PointSet<dim,Scalar>::end() const
{
  return points.end();
}
}
