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

size_t Graph::degree(const size_t v) const
{
  assert( v < nv() && "Invalid vertex index");
  return _adj[v].size();
}

Edge & Graph::edge(const size_t e)
{
  assert( e < ne() && "Invalid edge index" );
  return _edges[e];
}

}  // end namespace algorithms
