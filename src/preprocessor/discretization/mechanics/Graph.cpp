#include "Graph.hpp"
#include "angem/utils.hpp"
#include <cassert>

namespace algorithms {

Graph::Graph(const size_t nv) : _adj(nv)
{}

void Graph::add(const size_t v1, const size_t v2)
{
  size_t const e = ne();
  _edges.emplace_back(v1, v2);
  _adj[v1].push_back(e);
  _adj[v2].push_back(e);
}

std::vector<Edge const*> Graph::adj(const size_t v) const
{
  assert( v < nv() && "Invalid vertex index");
  std::vector<Edge const *> ans;
  ans.reserve(degree(v));
  for (size_t e : _adj[v])
    ans.push_back(&_edges[e]);
  return ans;
}

size_t Graph::degree(const size_t v) const
{
  assert( v < nv() && "Invalid vertex index");
  return _adj[v].size();
}

std::vector<size_t> const &  Graph::adj_idx(const size_t v) const
{
  assert( v < nv() && "Invalid vertex index");
  return _adj[v];
}

Edge & Graph::edge(const size_t e)
{
  assert( e < ne() && "Invalid edge index" );
  return _edges[e];
}

void Graph::reorder_vertex_edges(size_t v, std::vector<size_t> & order)
{
  angem::reorder<size_t, size_t>(_adj[v], order);
}


}  // end namespace algorithms
