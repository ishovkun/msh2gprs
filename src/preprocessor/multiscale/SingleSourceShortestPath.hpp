#pragma once
#include "EdgeWeightedDigraph.hpp"
#include <limits>  // numeric_limits
#include <algorithm>  // reverse
#include <iostream>

namespace algorithms {

class SingleSourceShortestPath {
 public:
  virtual ~SingleSourceShortestPath() = default;

  double distanceTo(const size_t vertex) const {return _dist_to[vertex];}

  std::vector<DirectedEdge> pathTo(const size_t vertex) const
  {
    std::vector<DirectedEdge> result;
    auto * p_edge = &_edge_to[vertex];
    for (; p_edge->from() != _source; p_edge = &_edge_to[p_edge->from()])
      result.push_back(*p_edge);
    result.push_back(*p_edge);  // edge to the source vertex
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
    for (const auto & edge : pathTo(vertex))
      std::cout << edge.from() << "->" << edge.to() << " (" << edge.weight() << ") ";
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
      _edge_to[edge.to()] = edge;
    }
  }

  const EdgeWeightedDigraph & _gr;
  const size_t _source;
  std::vector<DirectedEdge> _edge_to;
  std::vector<double> _dist_to;
};

}  // end namespace algorithms
