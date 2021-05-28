#pragma once

#include "DirectedEdge.hpp"
#include <list>
#include <vector>
#include <stdexcept>  // invalid_argument
#include <istream>    // stream

namespace algorithms {

class EdgeWeightedDigraph {
 public:
  EdgeWeightedDigraph(const size_t n_vertices) noexcept
      : _adj(n_vertices)
  {}

  EdgeWeightedDigraph(std::istream & in)
  {
    size_t nvert, nconn, v1, v2;
    double weight;
    in >> nvert;
    _adj.resize(nvert);
    in >> nconn;
    for (size_t i=0; i<nconn; ++i)
    {
      in >> v1;
      in >> v2;
      in >> weight;
      add(DirectedEdge(v1, v2, weight));
    }


  }

  void add(const DirectedEdge & edge)
  {
    checkEdgeValidity_(edge);
    _edges.push_back(edge);
    _adj[edge.from()].push_back(_edges.size() - 1);
  }

  std::vector<DirectedEdge const *> adj(const size_t v) const
  {
    checkVertexValidity_(v);
    std::vector<DirectedEdge const*> result;
    result.reserve(degree(v));
    for (size_t e : _adj[v])
      result.push_back( &_edges[e] );
    return result;
  }

  std::vector<DirectedEdge*> adj(const size_t v)
  {
    checkVertexValidity_(v);
    std::vector<DirectedEdge*> result;
    result.reserve(degree(v));
    for (size_t e : _adj[v])
      result.push_back( &_edges[e] );
    return result;
  }

  // number of outgoing edges
  inline size_t degree(size_t v) const noexcept {return _adj[v].size();}
  inline std::vector<DirectedEdge> const & edges() const noexcept {return _edges;}
  inline std::vector<DirectedEdge> & edges() noexcept {return _edges;}

  size_t n_vertices() const noexcept {return _adj.size();}
  size_t n_edges() const noexcept {return _edges.size();}

  virtual ~EdgeWeightedDigraph() = default;

 protected:
  void checkVertexValidity_(const size_t v) const
  {
    if (v >= n_vertices())
      throw std::invalid_argument("vertex index > n_vertices");
  }

  void checkEdgeValidity_(const DirectedEdge & e) const
  {
    checkVertexValidity_(e.from());
    checkVertexValidity_(e.to());
  }


  std::vector<DirectedEdge> _edges;    // storage for edges
  std::vector<std::list<size_t>> _adj;  // vertex -> list of adjacent edge indices
};

}  // end namespace algorithms
