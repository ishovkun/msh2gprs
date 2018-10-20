#pragma once

#include <Point.hpp>
#include <Polygon.hpp>
#include <PointSet.hpp>
// #include <PolyGroup.hpp>


namespace angem
{

template <typename Scalar>
class SurfaceMesh // : PolyGroup<Scalar>
{
 public:
  SurfaceMesh(const double tol = 1e-6,
              const std::size_t max_edges = 1e8);
  void insert(const Polygon<Scalar> & poly);
  bool empty() const;

  // merge polygon with its larger neighbor
  void merge_element(const std::size_t element);

  PointSet<3,Scalar>                    vertices;
  std::vector<std::vector<std::size_t>> polygons;  // indices
  // std::vector<std::size_t> neighbors;

  // hash of two vert indices -> vector polygons
  // essentially edge -> neighbors
  std::unordered_map<std::size_t, std::vector<std::size_t>> map_edge_neighbors;

 private:
  std::size_t hash_value(const std::size_t ind1,
                         const std::size_t ind2) const
  {
    if (ind1 < ind2)
      return max_edges*ind1 + ind2;
    else
      return max_edges*ind2 + ind1;
  }


  Scalar tol;
  std::size_t max_edges;  // used to hash a int-int pair
};


template <typename Scalar>
SurfaceMesh<Scalar>::SurfaceMesh(const double tol,
                                 const std::size_t max_edges)
    :
    tol(tol),
    max_edges(max_edges)
{
  assert(tol > 0);
}


template <typename Scalar>
void SurfaceMesh<Scalar>::insert(const Polygon<Scalar> & poly)
{
  // this part is a copy of PolyGroup add method
  // append vertices, compute vertex indices, and append these indices
  std::vector<std::size_t> indices;
  const std::vector<Point<3,Scalar>> & points = poly.get_points();
  for (const auto & p : points)
  {
    const std::size_t ind = vertices.insert(p);
    indices.push_back(ind);
  }
  polygons.push_back(indices);

  // now we need to compute edges
  for (std::size_t i=0; i<points.size(); ++i)
  {
    std::vector<std::size_t> * p_edge;
    std::size_t i1, i2;
    if (i < points.size() - 1)
    {
      i1 = vertices.find(points[i]);
      i2 = vertices.find(points[i+1]);
    }
    else
    {
      i1 = vertices.find(points[i]);
      i2 = vertices.find(points[0]);
    }

    if (i1 == i2)
      throw std::runtime_error("what kind of polygon is that?");

    std::size_t hash = hash_value(i1, i2);
    auto iter = map_edge_neighbors.find(hash);
    if (iter != map_edge_neighbors.end())
      (iter->second).push_back(polygons.size());
    else
      map_edge_neighbors.insert({ {hash, {polygons.size()}} });
  }

}


template <typename Scalar>
void SurfaceMesh<Scalar>::merge_element( const std::size_t element )
{
  // make iterator for faces, i'm sick of this code
  // for (std::size_t i=0; i<points.size(); ++i)
  //   if (i < points.size() - 1)
  //   else

  Polygon<Scalar> poly(vertices, polygons(element));
  std::vector<Point<3,Scalar>> &points = poly.get_points();
  // find first neighbor
  for (std::size_t i=0; i<points.size(); ++i)
  {
    std::size_t i1, i2;
    if (i < points.size() - 1)
    {
      i1 = vertices.find(points[i]);
      i2 = vertices.find(points[i+1]);
    }
    else
    {
      i1 = vertices.find(points[i]);
      i2 = vertices.find(points[0]);
    }

    std::size_t hash = hash_value(i1, i2);
    auto iter = map_edge_neighbors.find(hash);
    if (iter != map_edge_neighbors.end())  // found!!!
    {

    }
    else
      continue;

  }
}

}
