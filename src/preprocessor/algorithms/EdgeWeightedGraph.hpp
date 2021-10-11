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
  // check if edge between vertices u and v exists
  bool has_edge(size_t u, size_t v) const;
  // edges incident to i
  std::vector<UndirectedEdge const*> adj(size_t i) const;
  // non-const edges incident to i
  std::vector<UndirectedEdge*> adj(size_t i);
  // all edges in this graph
  const std::vector<UndirectedEdge> & edges() const {return _edges;}
  // all edges in this graph
  std::vector<UndirectedEdge> & edges() {return _edges;}
  // number of vertices
  size_t n_vertices() const {return _adj.size();}
  // number of edges
  size_t n_edges() const {return _edges.size();}
  // prints the graph
  friend  std::ostream & operator<<(std::ostream & os, const EdgeWeightedGraph & gr);

  std::vector<UndirectedEdge> _edges;
  std::vector<std::list<size_t>> _adj;
};


}  // end namespace algorithms
