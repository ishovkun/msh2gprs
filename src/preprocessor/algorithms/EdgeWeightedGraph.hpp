#pragma once

#include "UndirectedEdge.hpp"
#include <vector>
#include <list>

namespace algorithms {

class EdgeWeightedGraph
{
 public:
  // constructor
  EdgeWeightedGraph(const size_t nv);
  // add weighted edge e to this graph
  void add(const UndirectedEdge & e);
  // edges incident to i
  std::vector<UndirectedEdge const*> adj(size_t i) const;
  // non-const edges incident to i
  std::vector<UndirectedEdge*> adj(size_t i);
  // all edges in this graph
  const std::vector<UndirectedEdge> & edges() const {return _edges;}
  // number of vertices
  size_t n_vertices() const {return _adj.size();}
  // number of edges
  size_t n_edges() const {return _edges.size();}
  //
  friend  std::ostream & operator<<(std::ostream & os, const EdgeWeightedGraph & gr);

  std::vector<UndirectedEdge> _edges;
  std::vector<std::list<size_t>> _adj;
};


}  // end namespace algorithms
