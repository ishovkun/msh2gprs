#pragma once
#include "FlowEdge.hpp"
#include <vector>
#include <algorithm>  // transform
#include <istream>    // istream
#include <iostream>   // debug


namespace algorithms {

class FlowNetwork {
 public:
  FlowNetwork(size_t nv)
      : _adj(nv, std::vector<size_t>())
  {}

  FlowNetwork(std::istream & in)
  {
    size_t nvert, nconn, v1, v2;
    double weight;
    in >> nvert;

    if (nvert == 0)
      throw std::invalid_argument("FN has zero vertices");
    if (nvert > 1e6)
      throw std::invalid_argument("FN has too many vertices");

    _adj.resize(nvert);
    in >> nconn;
    for (size_t i=0; i<nconn; ++i)
    {
      in >> v1;
      in >> v2;
      in >> weight;
      validateVertex_(v1);
      validateVertex_(v2);
      addEdge(FlowEdge(v1, v2, weight));
    }
  }

  virtual ~FlowNetwork() = default;

  void addEdge(const FlowEdge & edge)
  {
    validateVertex_(edge.from());
    validateVertex_(edge.to());
    _edges.push_back(edge);
    const size_t iedge = n_edges() - 1;
    _adj[edge.from()].push_back(iedge);
    _adj[edge.to()].push_back(iedge);
  }

  std::vector<const FlowEdge*> adj(const size_t v) const
  {
    std::vector<const FlowEdge*> result(_adj[v].size(), nullptr);
    std::transform(_adj[v].begin(), _adj[v].end(), result.begin(),
                   [this](const size_t iedge) {return &_edges[iedge];});
    return result;
  }

  const std::vector<FlowEdge> & edges() const noexcept {return _edges;}
  size_t n_vertices() const noexcept {return _adj.size();}
  size_t n_edges() const noexcept {return _edges.size();}

 protected:
  inline void validateVertex_(const size_t v) {assert( v < n_vertices() && "Invalid vertex" );}

  std::vector<std::vector<size_t>> _adj;
  std::vector<FlowEdge> _edges;
};

}  // end namespace algorithms
