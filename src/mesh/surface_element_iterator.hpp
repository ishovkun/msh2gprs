#pragma once

#include "angem/PointSet.hpp"
#include <edge_iterator.hpp>
#include <surface_mesh_methods.hpp>

namespace mesh
{

template <typename Scalar>
class surface_element_iterator
{
 public:
  surface_element_iterator(const std::size_t                       ielement,
                           angem::PointSet<3,Scalar>             & vertices,
                           std::vector<std::vector<std::size_t>> & polygons,
                           EdgeMap                               & map_edges,
                           std::vector<int>                & grid_markers);

  // assignment operator
  surface_element_iterator<Scalar> & operator=(const surface_element_iterator<Scalar> & other);
  // comparison operators
  bool operator==(const surface_element_iterator<Scalar> & other) const;
  bool operator!=(const surface_element_iterator<Scalar> & other) const;
  // incrementing/decrementing
  surface_element_iterator<Scalar> & operator++();
  surface_element_iterator<Scalar> & operator--();
  // Getters
  inline std::size_t index() const {return ielement;}
  inline int marker() const { return _grid_markers[ielement]; }
  std::vector<surface_element_iterator<Scalar>> neighbors() const;
  std::vector<edge_iterator<Scalar>> edges() const;
  const std::vector<std::size_t> & vertices() const { return polygons[ielement]; }
  std::vector<std::size_t> & vertices() { return polygons[ielement]; };
  std::vector<angem::Point<3,Scalar>> vertex_coordinates() const;

 private:
  std::size_t                             ielement;
  angem::PointSet<3,Scalar>             & mesh_vertices;
  std::vector<std::vector<std::size_t>> & polygons;
  EdgeMap                               & map_edges;
  std::vector<int>                & _grid_markers;
};


template <typename Scalar>
surface_element_iterator<Scalar>::
surface_element_iterator(const std::size_t                       ielement,
                         angem::PointSet<3,Scalar>             & vertices,
                         std::vector<std::vector<std::size_t>> & polygons,
                         EdgeMap                               & map_edges,
                         std::vector<int>                      & grid_markers)
    : ielement(ielement),
      mesh_vertices(vertices),
      polygons(polygons),
      map_edges(map_edges),
      _grid_markers(grid_markers)
{}


template <typename Scalar>
surface_element_iterator<Scalar> & surface_element_iterator<Scalar>::operator++()
{
  ielement++;
  return *this;
}


template <typename Scalar>
surface_element_iterator<Scalar> & surface_element_iterator<Scalar>::operator--()
{
  ielement--;
  return *this;
}


template <typename Scalar>
bool surface_element_iterator<Scalar>::
operator==(const surface_element_iterator<Scalar> & other) const
{
  return ielement == other.ielement;
}


template <typename Scalar>
bool surface_element_iterator<Scalar>::
operator!=(const surface_element_iterator<Scalar> & other) const
{
  return !(*this == other);
}


template <typename Scalar>
std::vector<edge_iterator<Scalar>>
surface_element_iterator<Scalar>::edges() const
{
  std::vector<edge_iterator<Scalar>> edges;

  std::vector<angem::Point<3,Scalar>> points;
  for (const auto & vertex : polygons[ielement])
    points.push_back(mesh_vertices.points[vertex]);

  for (std::size_t i=0; i<points.size(); ++i)
  {
    std::size_t i1, i2;
    if (i < points.size() - 1)
    {
      i1 = mesh_vertices.find(points[i]);
      i2 = mesh_vertices.find(points[i+1]);
    }
    else
    {
      i1 = mesh_vertices.find(points[i]);
      i2 = mesh_vertices.find(points[0]);
    }

    const auto hash = hash_value({i1, i2});
    edges.push_back(edge_iterator<Scalar>(map_edges.find(hash), mesh_vertices));
  }
  return edges;
}


template <typename Scalar>
std::vector<angem::Point<3,Scalar>> surface_element_iterator<Scalar>::vertex_coordinates() const
{
  std::vector<angem::Point<3,Scalar>> result;
  result.reserve(polygons[ielement].size());

  for (std::size_t ivertex : polygons[ielement])
    result.push_back(mesh_vertices[ivertex]);

  return result;
}


}  // end namespace
