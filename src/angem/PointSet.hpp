#pragma once
#include <Point.hpp>

#include <unordered_map>
#include <algorithm>

// 128 ints for point hashing
#ifdef USE_BOOST
#include <boost/multiprecision/cpp_int.hpp>
#else
#include <uint256/uint128_t.h>
#endif

namespace angem
{

#ifdef USE_BOOST
using pset_hash_type = boost::multiprecision::uint128_t;
#else
using pset_hash_type = uint128_t;
#endif


/* Set of points
 * Lookup and insertion in O(1) time due to hashing
 * the idea of hashing is to divide 3D space into cartesian blocks
 * id_i = (point[i] - minimum_number)/tolerance
 * hash = id_x + id_y*nx + id_z*(nx*ny)
 * Note: I use 128-bit ints and doubles for this calculations
 * might be redundant
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
  std::size_t find(const Point<dim,Scalar> &p) const;

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
  bool empty() const {return points.empty();}

  // variables
  std::vector<Point<dim,Scalar>> points;
 private:
  pset_hash_type hash_value(const Point<dim,Scalar> &) const;

  // variables
  // bounds
  Point<dim,long double> lower, upper;
  std::unordered_map<pset_hash_type,std::size_t> index_map;
  // search tolerance
  double tol;
};


template<int dim, typename Scalar>
std::size_t PointSet<dim,Scalar>::find(const Point<dim,Scalar> &p) const
{
  const pset_hash_type ha = hash_value(p);
  const std::unordered_map<pset_hash_type,std::size_t>::const_iterator it = index_map.find(ha);

  if (it == index_map.end())
    return size();
  else
    return it->second;
}


template<int dim, typename Scalar>
std::size_t PointSet<dim,Scalar>::insert(const Point<dim,Scalar> &p)
{
  const pset_hash_type ha = hash_value(p);
  const std::unordered_map<pset_hash_type,std::size_t>::const_iterator it = index_map.find(ha);
  if (it == index_map.end())
  {
    const std::size_t ind = size();
    index_map.insert({ha, ind});
    points.push_back(p);
    return ind;
  }
  else
    return it->second;
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


// helper function for pset_hash_type output (debugging)
// std::ostream&
// operator<<( std::ostream& dest, pset_hash_type_t value )
// {
//   std::ostream::sentry s( dest );
//   if ( s ) {
//     __uint128_t tmp = value < 0 ? -value : value;
//     char buffer[ 128 ];
//     char* d = std::end( buffer );
//     do
//     {
//       -- d;
//       *d = "0123456789"[ tmp % 10 ];
//       tmp /= 10;
//     } while ( tmp != 0 );
//     if ( value < 0 ) {
//       -- d;
//       *d = '-';
//     }
//     int len = std::end( buffer ) - d;
//     if ( dest.rdbuf()->sputn( d, len ) != len ) {
//       dest.setstate( std::ios_base::badbit );
//     }
//   }
//   return dest;
// }


template<int dim, typename Scalar>
pset_hash_type PointSet<dim,Scalar>::hash_value(const Point<dim,Scalar> &p) const
{
  assert((upper-lower).norm() > tol);

  pset_hash_type result = 0;
  pset_hash_type mult = 1;

  for (short i=0; i<dim; ++i)
  {
    const long double ixd =
        (static_cast<long double>(p[i]) -
        static_cast<long double>(lower[i])) /
        static_cast<long double>(tol) ;

    if (ixd < 0 or ixd > (upper[i] - lower[i]) / tol)
    {
      std::cout.precision(12);
      // std::cout << "ixd = " << ixd << std::endl;
      std::cout << "weird point = " << p << std::endl;
      throw std::out_of_range("Point cannot be added to set");
    }

    const pset_hash_type ix = static_cast<pset_hash_type> (std::round(ixd));

    result += ix * mult;

    const pset_hash_type nx = static_cast<pset_hash_type>((upper[i] - lower[i]) / tol);
    mult *= nx;
  }

  return result;
}

}  // end namespace
