#pragma once

#include <angem/Point.hpp>
#include <angem/PointSet.hpp>
#include <edge_iterator.hpp>

#include <unordered_map>
#include <vector>


namespace mesh
{
using EdgeMap = std::unordered_map<std::size_t, std::vector<std::size_t>>;

template <typename Scalar>
class const_edge_iterator
{
 public:
  const_edge_iterator(EdgeMap::const_iterator         it,
                      const angem::PointSet<3,Scalar> & vertices);
  // comparison
  bool operator==(const edge_iterator<Scalar> & other) const;
  bool operator!=(const edge_iterator<Scalar> & other) const;
  bool operator==(const const_edge_iterator<Scalar> & other) const;
  bool operator!=(const const_edge_iterator<Scalar> & other) const;
  // incrementing
  const_edge_iterator<Scalar> & operator++();

  // Getters
  inline const std::vector<std::size_t> & neighbors() const {return map_it->second;}
  std::pair<std::size_t, std::size_t> vertex_indices() const;

 private:
  EdgeMap::const_iterator map_it;
  const angem::PointSet<3,double> & mesh_vertices;
};


template <typename Scalar>
const_edge_iterator<Scalar>::
const_edge_iterator(EdgeMap::const_iterator           it,
                    const angem::PointSet<3,Scalar> & vertices)
    :
    map_it(it),
    mesh_vertices(vertices)
{}


template <typename Scalar>
const_edge_iterator<Scalar> & const_edge_iterator<Scalar>::operator++()
{
  map_it++;
  return *this;
}


template <typename Scalar>
bool const_edge_iterator<Scalar>::
operator==(const edge_iterator<Scalar> & other) const
{
  if (map_it != other.map_it)
    return false;
  else
    return true;
}


template <typename Scalar>
bool const_edge_iterator<Scalar>::
operator==(const const_edge_iterator<Scalar> & other) const
{
  if (map_it != other.map_it)
    return false;
  else
    return true;
}


template <typename Scalar>
bool const_edge_iterator<Scalar>::
operator!=(const edge_iterator<Scalar> & other) const
{
  return !(*this == other);
}


template <typename Scalar>
bool const_edge_iterator<Scalar>::
operator!=(const const_edge_iterator<Scalar> & other) const
{
  return !(*this == other);
}


template <typename Scalar>
std::pair<std::size_t, std::size_t>
const_edge_iterator<Scalar>::vertex_indices() const
{
  return invert_hash(map_it->first);
}


}
