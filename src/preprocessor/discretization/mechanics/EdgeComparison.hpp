#pragma once

#include "ConnectionMap.hpp"
#include <vector>
#include <list>
#include <algorithm> // std min

namespace discretization {

namespace edgecmp {

class Edge
{
 public:
  Edge() : _first(0), _second(0){}
  Edge(const size_t first, const size_t second)
      : _first(first), _second(second) {}

  Edge(const Edge & e)
  {
    _first = e.either();
    _second = e.other(_first);
  }

  Edge & operator=(const Edge & e)
  {
    _first = e.either();
    _second = e.other(_first);
    return *this;
  }

  size_t either() const { return _first; }
  size_t other(const size_t v) const {return (v == _first) ? _second : _first;}

  bool operator==(const Edge & e) const
  {
    if ( min_() == e.min_() && max_() == e.max_() )
      return true;
    else return false;
  }

 private:
  size_t min_() const { return std::min(_first, _second); }
  size_t max_() const { return std::max(_first, _second); }

  size_t _first, _second;
};

class EdgeComparison
{
 public:

  static hash_algorithms::ConnectionMap<std::vector<Edge>>
  get_edges(const std::vector<std::list<size_t>> & vertex_face_markers)
  {
    hash_algorithms::ConnectionMap<std::vector<Edge>> cm;
    for (size_t v=0; v<vertex_face_markers.size(); ++v)
      for (size_t w=v+1; w<vertex_face_markers.size(); ++w)
      {
        // if face pair matches then there is an edge between verttices
        for (const Edge & ev : get_permutations(vertex_face_markers[v]))
          for (const Edge & ew : get_permutations(vertex_face_markers[w]))
            if ( ev == ew )
            {
              const size_t m1 = ev.either();
              const size_t m2 = ev.other(m1);
              if (!cm.contains(m1, m2))
              {
                const size_t data_index = cm.insert( m1, m2 );
                cm.get_data(data_index).push_back(Edge(v , w));
              }
              else
              {
                cm.get_data(m1, m2).push_back(Edge(v , w));
              }
            }
      }

    return cm;
  }

 private:
  static std::vector<Edge> get_permutations(const std::list<size_t> & values)
  {
    std::vector<Edge> result;
    for (auto iit = values.begin(); iit != values.end(); ++iit)
    {
      auto jit = iit; ++jit;
      for (; jit != values.end(); ++jit)
        result.push_back(Edge( *iit, *jit ));
    }
    return result;
  }
};


}  // end namespace edgecmp



}  // end namespace discretization
