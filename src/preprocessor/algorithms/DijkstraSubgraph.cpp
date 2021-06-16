#include "DijkstraSubgraph.hpp"
#include <limits>
#include <cassert>

namespace algorithms {

DijkstraSubgraph::DijkstraSubgraph(size_t source, EdgeWeightedGraph const & g, std::vector<bool> const &bitset)
    : _source(source)
    , _g(g)
    , _bitset(bitset)
{
  assert( _bitset[_source] && "source must be in subgraph" );

  for (size_t i = 0; i < _bitset.size(); ++i)
    if (_bitset[i])
      _path.insert( { i, PathVertex{nullptr, std::numeric_limits<double>::max()}  });

  _path[_source].dist = 0;
  _pq.push( { _source, 0.0 } );

  while (!_pq.empty()) {
    const size_t u = _pq.top().u;
    _pq.pop();
    for (auto const e : _g.adj(u))
      if (_bitset[e->other(u)]) {
        relax_(*e, u);
      }
  }
}

void DijkstraSubgraph::relax_(UndirectedEdge const & e, size_t from)
{
  size_t const u = from, v = e.other(u);
  double const new_dist = _path[u].dist + e.weight();
  if ( new_dist < _path[v].dist ) {
    _path[v] = { &e, new_dist };
    _pq.push({ v, new_dist });
  }
}

double DijkstraSubgraph::distanceTo(size_t u) const
{
  assert( _bitset[u] && "vertex is not in subgraph" );
  return _path.find(u)->second.dist;
}

std::vector<UndirectedEdge const*> DijkstraSubgraph::pathTo(size_t u) const
{
  assert( _bitset[u] && "vertex is not in subgraph" );
  std::vector<UndirectedEdge const*> ans;
  if (u == _source) return ans;

  auto * e = _path.find(u)->second.edge;
  assert( e && " path to this vertex apparently does not exist" );

  for (; e->other(u) != _source; e = _path.find( e->other(u) )->second.edge ) {
    ans.push_back(e);
    u = e->other(u);
  }
  return ans;
}

}  // end namespace algorithms
