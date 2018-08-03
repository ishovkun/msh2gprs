#pragma once

#include "Point.hpp"

namespace angem
{

template <int dim, typename Scalar>
class Basis
{
 public:
  Basis();
  Basis(const std::vector<Point<dim,Scalar>> & vecs);

  // setters
  void set_data(const std::vector<Point<dim,Scalar>> & vecs);

  Point<dim,Scalar> & operator[] (int i);
  const Point<dim,Scalar> & operator() (int i) const;

  Point<dim,Scalar> transform(const Point<dim,Scalar> & p) const;
  bool is_empty() const;

 private:
  std::vector<Point<dim,Scalar>> vectors;
};


template <int dim, typename Scalar>
Basis<dim,Scalar>::Basis()
{
  vectors.resize(dim);
}


template <int dim, typename Scalar>
Basis<dim,Scalar>::Basis(const std::vector<Point<dim,Scalar>> & vecs)
    :
    vectors(vecs)
{
  assert(vecs.size() == dim);
}


template <int dim, typename Scalar>
void
Basis<dim,Scalar>::set_data(const std::vector<Point<dim,Scalar>> & vecs)
{
  assert(vecs.size() == dim);
  vectors = vecs;
}


template <int dim, typename Scalar>
Point<dim,Scalar> &
Basis<dim,Scalar>::operator[](int i)
{
  assert(i < dim);
  return vectors[i];
}


template <int dim, typename Scalar>
const Point<dim,Scalar> &
Basis<dim,Scalar>::operator()(int i) const
{
  assert(i < dim);
  return vectors[i];
}


template <int dim, typename Scalar>
Point<dim,Scalar>
Basis<dim,Scalar>::transform(const Point<dim,Scalar> & p) const
{
  assert(!is_empty());
  Point<dim,Scalar> result;
  for (int i=0; i<dim; ++i)
    result[i] = p.dot(vectors[i]);
  return result;
}


template <int dim, typename Scalar>
bool
Basis<dim,Scalar>::is_empty() const
{
  Scalar sum = static_cast<Scalar>(0);
  for (int i=0; i<dim; ++i)
    sum += vectors[i].norm();
  if (fabs(sum - static_cast<Scalar>(3)) < 1e-16)
    return true;
  else
    return false;

}
}  // end namespace
