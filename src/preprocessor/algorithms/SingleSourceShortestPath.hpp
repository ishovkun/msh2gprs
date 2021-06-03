#pragma once
#include "EdgeWeightedDigraph.hpp"
#include <limits>  // numeric_limits
#include <algorithm>  // reverse
#include <iostream>
#include <cassert>

namespace algorithms {

class SingleSourceShortestPath {
 public:
  virtual ~SingleSourceShortestPath() = default;

  size_t source() const noexcept {return _source;}

  double distanceTo(const size_t vertex) const {return _dist_to[vertex];}

  std::vector<DirectedEdge const*> pathTo(const size_t vertex) const
  {
    std::vector<DirectedEdge const*> result;
    if (vertex == _source) return result;

    assert( vertex <  _gr.n_vertices() && "invalid vertex requested");

    auto * edge = _edge_to[vertex];

    assert( edge && " path to this vertex apparently does not exist" );

    for (; edge->from() != _source; edge = _edge_to[edge->from()])
      result.push_back(edge);
    result.push_back(edge);  // edge to the source vertex
    std::reverse(result.begin(), result.end());
    return result;
  }

  void printPathTo(const size_t vertex) const
  {
    if (vertex == _source)
    {
      std::cout << 0 << std::endl;
      return;
    }
    for (const auto * edge : pathTo(vertex))
      std::cout << edge->from() << "->" << edge->to() << " (" << edge->weight() << ") ";
    std::cout << std::endl;
  }

 protected:
  SingleSourceShortestPath(const EdgeWeightedDigraph & gr, const size_t source)
      : _gr(gr), _source(source)
  {}

  virtual void relax_(const DirectedEdge & edge)
  {
    const double new_dist = _dist_to[edge.from()] + edge.weight();
    if ( new_dist < _dist_to[edge.to()] )
    {
      _dist_to[edge.to()] = new_dist;
      _edge_to[edge.to()] = &edge;
    }
  }

  EdgeWeightedDigraph const & _gr;
  const size_t _source;
  std::vector<DirectedEdge const *> _edge_to;
  std::vector<double> _dist_to;
};

}  // end namespace algorithms
