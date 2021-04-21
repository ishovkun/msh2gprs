#pragma once

#include <vector>
#include <fstream>
#include <iostream>
#include <stack>
#include <limits>

namespace algorithms {

using std::iostream;
using std::fstream;

class Edge
{
 public:
  explicit Edge(size_t v, size_t w) : v1(v), v2(w) {}
  size_t either() const {return v1;}
  size_t other(size_t v) const {return (v == v1) ? v2 : v1;}

 private:
  size_t v1, v2;
};

class Graph
{
 public:
  Graph(const size_t nv);
  void add(const size_t v1,  const size_t v2);
  std::vector<Edge const*> adj(const size_t v) const;
  std::vector<size_t> const & adj_idx(const size_t) const;
  Edge & edge(const size_t idx);

  size_t nv() const {return _adj.size();}
  size_t ne() const {return _edges.size();}
  size_t degree(const size_t v) const;

  friend std::ostream & operator<<(std::ostream & os, const Graph & gr);

 protected:
  std::vector<std::vector<size_t>> _adj;
  std::vector<Edge> _edges;
};


}  // end namespace algorithms
