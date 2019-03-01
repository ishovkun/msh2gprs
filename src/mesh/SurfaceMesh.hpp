#pragma once

#include <angem/Point.hpp>
#include <angem/Polygon.hpp>
#include <angem/PointSet.hpp>
#include <angem/Exceptions.hpp>
#include <surface_mesh_methods.hpp>
#include <edge_iterator.hpp>
#include <const_edge_iterator.hpp>
#include <surface_element_iterator.hpp>



namespace mesh
{

using namespace angem;

template <typename Scalar>
class SurfaceMesh // : PolyGroup<Scalar>
{
 public:
  // SurfaceMesh(const double tol = 1e-6, const std::size_t max_edges = 1e8);
  SurfaceMesh(const double tol = 1e-6);
  // add new polygon (element) to the mesh. no check for duplicates
  // returns new element index
  std::size_t insert(const angem::Polygon<Scalar> & poly);
  // true if vector of polygons is empty
  bool empty() const {return polygons.empty();}

  // GETTERS
  // get number of elements
  std::size_t n_polygons() const {return polygons.size();}
  // get vector of neighbor indices
  std::vector<std::size_t> get_neighbors( const std::size_t ielement ) const;
  // vector of indices of cells neighboring an edge
  std::vector<std::size_t> get_neighbors( const Edge & edge ) const;
  // get vector of index pairs representing edges
  std::vector<Edge> get_edges( const std::size_t ielement ) const;
  // delete an element from the mesh
  void delete_element(const std::size_t ielement);
  // merge polygon with its largest neighbor
  // complexity linear to number of faces (edges)
  // does not remove vertices when no connected neighbors
  std::size_t merge_element(const std::size_t element);
  Scalar minimum_edge_size() const;

  // ITERATORS
  edge_iterator<Scalar> create_edge_iterator(EdgeMap::iterator it)
  {
    return edge_iterator<Scalar>(it, vertices);
  }
  edge_iterator<Scalar> begin_edges() {return create_edge_iterator(map_edges.begin());}
  edge_iterator<Scalar> end_edges() {return create_edge_iterator(map_edges.end());}

  const_edge_iterator<Scalar> create_edge_iterator(EdgeMap::const_iterator it) const
  {
    return const_edge_iterator<Scalar>(it, vertices);
  }
  // create edge iterator starting from first edge in hash_map
  const_edge_iterator<Scalar> begin_edges() const {return create_edge_iterator(map_edges.begin());}
  // create end iterator; @TODO: should be const iterator?
  const_edge_iterator<Scalar> end_edges() const {return create_edge_iterator(map_edges.end());}

  surface_element_iterator<Scalar> create_poly_iterator(const std::size_t ielement)
  {
    return surface_element_iterator<Scalar> (ielement, vertices,
                                             polygons, map_edges);
  }

  surface_element_iterator<Scalar> begin_polygons() {return create_poly_iterator(0);}
  surface_element_iterator<Scalar> end_polygons() {return create_poly_iterator(polygons.size());}


  PointSet<3,Scalar>                    vertices;
  std::vector<std::vector<std::size_t>> polygons;  // indices

  // hash of two vert indices -> vector polygons
  // essentially edge -> neighbor elements
  std::unordered_map<std::size_t, std::vector<std::size_t>> map_edges;

 private:
  // merge two elements haaving a comon edge
  // the edge is vertex indices
  // private cause no check if they are neighbors
  // jelement is the element to be merged into ielement
  // does not remove vertices when no connected neighbors
  void merge_elements(const std::size_t ielement,
                      const std::size_t jelement,
                      const Edge      & edge);
  // delete and replace neighbors with the new element -- useful for merging
  void delete_replace_connections(const std::size_t deleted_element,
                                  const std::size_t replacement_element);


  // used for joining two polygons
  // returns vector of contingent edges without the removed edge
  // the result is a continuous curve
  std::vector<Edge> remove_edge(const std::size_t ielement,
                                const Edge      & edge) const;
  // get vertex indices from vector of contingent edges
  static std::vector<std::size_t> get_vertices(const std::vector<Edge> & edges);


  Scalar tol;
};


template <typename Scalar>
// SurfaceMesh<Scalar>::SurfaceMesh(const double tol,
//                                  const std::size_t max_edges)
SurfaceMesh<Scalar>::SurfaceMesh(const double tol)
    :
    tol(tol)
{
  assert(tol > 0);
}


template <typename Scalar>
std::size_t SurfaceMesh<Scalar>::insert(const angem::Polygon<Scalar> & poly)
{
  // this part is a copy of PolyGroup add method
  // append vertices, compute vertex indices, and append these indices
  std::vector<std::size_t> indices;
  const std::vector<angem::Point<3,Scalar>> & points = poly.get_points();
  for (const auto & p : points)
  {
    const std::size_t ind = vertices.insert(p);
    indices.push_back(ind);
  }

  const std::size_t new_element_index = polygons.size();
  polygons.push_back(indices);

  // now we need to compute edges
  for (const auto & edge : poly.get_edges())
  {
    const std::size_t i1 = vertices.find(points[edge.first]);
    const std::size_t i2 = vertices.find(points[edge.second]);

    assert(i1 < vertices.size());
    assert(i2 < vertices.size());

    if (i1 == i2)
      throw std::runtime_error("what kind of polygon is that?");

    const std::size_t hash = hash_value(i1, i2);
    auto iter = map_edges.find(hash);
    if (iter != map_edges.end())
      (iter->second).push_back(new_element_index);
    else
      map_edges.insert({ {hash, {new_element_index}} });
  }
  return new_element_index;

}


template <typename Scalar>
std::vector<std::size_t>
SurfaceMesh<Scalar>::get_neighbors( const std::size_t ielement ) const
{
  std::vector<std::size_t> v_neighbors;
  std::vector<Edge> edges = get_edges(ielement);
  for (const Edge & edge : edges)
  {
    for (const auto & ineighbor : get_neighbors(edge))
    {
      if (ineighbor == ielement)
        continue;
      else
        v_neighbors.push_back(ineighbor);
    }
  }
  return v_neighbors;
}


template <typename Scalar>
std::size_t SurfaceMesh<Scalar>::merge_element( const std::size_t ielement )
{
  const std::vector<Edge> edges = get_edges(ielement);

  std::size_t max_iedge = 0;
  std::size_t max_neighbor = 0;
  std::size_t counter = 0;
  Scalar max_area = 0;
  for (const auto & edge : edges)
  {
    const std::vector<std::size_t> edge_neighbors = get_neighbors(edge);
    for (std::size_t element : edge_neighbors)
      if (element != ielement)
      {
        const angem::Polygon<Scalar> poly(vertices.points, polygons[element]);
        const Scalar area = poly.area();
        if (area > max_area)
        {
          max_area = area;
          max_neighbor = element;
          max_iedge = counter;
        }
      }
    counter++;
  }

  // @TODO: eliminate vertices if adjacent edges are parallel
  // Note: this should involve creating new edges
  merge_elements(max_neighbor, ielement, edges[max_iedge]);
  return max_neighbor;
}


template <typename Scalar>
void SurfaceMesh<Scalar>::delete_element(const std::size_t element)
{
  // delete polygon and reduce all polygons
  // in map_edges with higher indices by 1
  polygons.erase(polygons.begin() + element);
  for (auto iter = map_edges.begin();
       iter != map_edges.end(); ++iter)
  {
    for (std::size_t i=0; i<iter->second.size(); ++i)
    {
      std::size_t ielement = iter->second[i];
      if (ielement == element)
        iter->second.erase(iter->second.begin() + i);
      else if (ielement > element)
        iter->second[i]--;
    }
  }
}

template <typename Scalar>
void SurfaceMesh<Scalar>::delete_replace_connections(const std::size_t deleted_element,
                                                     const std::size_t replacement_element)
{
  assert(deleted_element != replacement_element);
  // delete polygon and reduce all polygons
  // in map_edges with higher indices by 1
  std::size_t repl = replacement_element;
  if (repl > deleted_element)
    repl--;

  polygons.erase(polygons.begin() + deleted_element);

  for (auto iter = map_edges.begin();
       iter != map_edges.end(); ++iter)
  {
    for (std::size_t i=0; i<iter->second.size(); ++i)
    {
      std::size_t & ielement = iter->second[i];
      if (ielement == deleted_element)
        ielement = repl;
      else if (ielement > deleted_element)
        ielement--;
    }
  }

}


template <typename Scalar>
std::vector<Edge> SurfaceMesh<Scalar>::remove_edge(const std::size_t element,
                                                   const Edge      & edge) const
{
  std::vector<Edge> result;
  // compose two vectors before and after the removed edge
  std::vector<Edge> edges_before, edges_after;
  bool before_switch = true;
  for (const auto & iedge : get_edges(element))
  {
    if ( hash_value(iedge) == hash_value(edge) )
      before_switch = false;
    else
    {
      if (before_switch)
        edges_before.push_back(iedge);
      else
        edges_after.push_back(iedge);
    }
  }
  // join two vectors
  for (const auto & iedge : edges_after)
    result.push_back(iedge);
  for (const auto & iedge : edges_before)
    result.push_back(iedge);

  return result;
}


// edges is an array of contingent edges
template <typename Scalar>
std::vector<std::size_t>
SurfaceMesh<Scalar>::get_vertices(const std::vector<Edge> & edges)
{
  std::vector<std::size_t> v_vertices;
  for (const Edge & edge : edges)
    v_vertices.push_back(edge.first);
  if (edges.back().second != edges.front().first)
    v_vertices.push_back(edges.back().second);
  return v_vertices;
}


template <typename Scalar>
void SurfaceMesh<Scalar>::merge_elements(const std::size_t ielement,
                                         const std::size_t jelement,
                                         const Edge      & edge)
{
  // connect edges in a ordered fashion and delete jelement
  // the resulting polygon might be non-convex
  // that's why we don't use default polygon renumbering
  // and do it manually
  std::vector<std::size_t> new_element;
  const std::vector<Edge> iedges = remove_edge(ielement, edge);
  const std::vector<Edge> jedges = remove_edge(jelement, edge);

  if (iedges.front().first == jedges.back().second)
  {
    for (const Edge & jedge : jedges)
      new_element.push_back(jedge.first);
    for (const Edge & iedge : iedges)
      new_element.push_back(iedge.first);
  }
  else if (iedges.front().first == jedges.front().first)
  {
    for (auto it = jedges.rbegin(); it != jedges.rend(); ++it)
      new_element.push_back( (*it).second );
    for (const Edge & iedge : iedges)
      new_element.push_back(iedge.first);
  }
  else
    throw NotImplemented("can this even happen?");

  // remove joint point if lies on a flat edge
  for (std::size_t i=0; i<new_element.size(); ++i)
    if (new_element[i] == edge.first or new_element[i] == edge.second)
    {
      angem::Point<3,Scalar> v1 = vertices[new_element[i+1]] - vertices[new_element[i]];
      angem::Point<3,Scalar> v2 = vertices[new_element[i]] - vertices[new_element[i-1]];
      if ( (v1.cross(v2)).norm() < 1e-10 )
      {
        new_element.erase(new_element.begin() + i);
      }
    }

  polygons[ielement] = new_element;
  delete_replace_connections(jelement, ielement);
}


template <typename Scalar>
std::vector<Edge>
SurfaceMesh<Scalar>::get_edges( const std::size_t ielement ) const
{
  std::vector<Edge> edges;

  std::vector<angem::Point<3,Scalar>> points;
  for (const auto & vertex : polygons[ielement])
    points.push_back(vertices.points[vertex]);

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

    edges.push_back({i1, i2});

  }
  return edges;
}


template <typename Scalar>
std::vector<std::size_t>
SurfaceMesh<Scalar>::get_neighbors( const Edge & edge ) const
{
  std::vector<std::size_t> v_neighbors;
  const std::size_t hash = hash_value(edge.first, edge.second);
  const auto iter = map_edges.find(hash);
  if (iter == map_edges.end())
    throw std::out_of_range("edge does not exist");

  return iter->second;
}


template <typename Scalar>
Scalar SurfaceMesh<Scalar>::minimum_edge_size() const
{
  Scalar min_size = std::numeric_limits<Scalar>::max();
  for (const auto & conn : map_edges)
  {
    const std::pair<std::size_t,std::size_t> verts = invert_hash(conn.first);
    const Scalar h = (vertices[verts.first] - vertices[verts.second]).norm();
    if (h < min_size)
      min_size = h;
  }

  return min_size;
}

}
