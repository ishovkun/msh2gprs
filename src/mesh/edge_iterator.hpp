#pragma once

#include "angem/Point.hpp"
#include "angem/PointSet.hpp"

#include <unordered_map>
#include <vector>


namespace mesh
{

using EdgeMap = std::unordered_map<std::size_t, std::vector<std::size_t>>;

template <typename Scalar>
class edge_iterator
{
 public:
  edge_iterator(EdgeMap::iterator          it,
                angem::PointSet<3,Scalar> & vertices);
  // comparison
  bool operator==(const edge_iterator<Scalar> & other) const;
  bool operator!=(const edge_iterator<Scalar> & other) const;
  // incrementing
  edge_iterator<Scalar> & operator++();

  // Getters
  // vector of surface mesh element indices
  inline const std::vector<std::size_t> & neighbors() const {return map_it->second;}
  std::pair<std::size_t, std::size_t> vertex_indices() const;
  std::pair<angem::Point<3,double>,angem::Point<3,double> > vertices() const;

 private:
  EdgeMap::iterator map_it;
  angem::PointSet<3,double> & mesh_vertices;
};


template <typename Scalar>
edge_iterator<Scalar>::edge_iterator(EdgeMap::iterator          it,
                                     angem::PointSet<3,Scalar> & vertices)
    :
    map_it(it),
    mesh_vertices(vertices)
{}


template <typename Scalar>
edge_iterator<Scalar> & edge_iterator<Scalar>::operator++()
{
  map_it++;
  return *this;
}


template <typename Scalar>
bool edge_iterator<Scalar>::operator==(const edge_iterator<Scalar> & other) const
{
  if (map_it != other.map_it)
    return false;
  else
    return true;
}


template <typename Scalar>
bool edge_iterator<Scalar>::operator!=(const edge_iterator<Scalar> & other) const
{
  return !(*this == other);
}


template <typename Scalar>
std::pair<std::size_t, std::size_t> edge_iterator<Scalar>::vertex_indices() const
{
  return invert_hash(map_it->first);
}


template <typename Scalar>
std::pair<angem::Point<3,double>,angem::Point<3,double> >
edge_iterator<Scalar>::vertices() const
{
  const auto ivertices = vertex_indices();
  return {mesh_vertices[ivertices.first], mesh_vertices[ivertices.second]};
}

}
