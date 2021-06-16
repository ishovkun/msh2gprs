#pragma once

#include "EdgeWeightedGraph.hpp"
#include <unordered_map>
#include <queue>

namespace algorithms {

/*
** Implements Dijkstra's shortes path in subgraph.
** Rutime complexity: E log(E)
*/
class DijkstraSubgraph {
 public:
  // bitset indicates subgraph vertices
  DijkstraSubgraph(size_t source, EdgeWeightedGraph const & g, std::vector<bool> const &bitset);
  double distanceTo(size_t u) const;
  std::vector<UndirectedEdge const*> pathTo(size_t u) const;

  virtual ~DijkstraSubgraph() = default;


 private:
  void relax_(UndirectedEdge const & e, size_t from);

  size_t const _source;
  EdgeWeightedGraph const & _g;
  std::vector<bool> const & _bitset;

  struct PathVertex {
    UndirectedEdge const * edge;
    double dist;
  };
  struct PrioritizedVertex {
    size_t u;
    double priority;
    bool operator>(PrioritizedVertex const &other) const { return priority > other.priority; }
  };

  std::unordered_map<size_t,PathVertex> _path;
  std::priority_queue<PrioritizedVertex, std::vector<PrioritizedVertex>, std::greater<PrioritizedVertex>> _pq;
};



}  // end namespace algorithms
