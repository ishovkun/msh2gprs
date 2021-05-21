#pragma once

#include <string>
#include <iostream>
#include <cassert>

namespace algorithms {

class UndirectedEdge
{
public:
  /* create an edge */
  UndirectedEdge(const size_t v, const size_t w, const double weight)
      : _first(v), _second(w), _weight(weight) {}
  // copy operator
  void operator=(const UndirectedEdge & e)
  {
    const size_t v = e.either();
    _first = v;
    _second = e.other(v);
    _weight = e.weight();
    // return ne;
  }
  // either endpoint
  inline size_t either() const noexcept {return _first;}
  // the endpoint that is not v
  size_t other(const size_t v) const {assert(v == _first || v == _second); return (v == _first) ? _second : _first;}
  // comparison by weight
  inline bool operator==(const UndirectedEdge & other) const noexcept {return weight() == other.weight();}
  // comparison by weight
  inline bool operator<(const UndirectedEdge & other) const noexcept {return weight() <  other.weight(); }
  // comparison by weight
  inline bool operator>(const UndirectedEdge & other) const noexcept {return weight() >  other.weight(); }
  // get the weight
  double weight() const {return _weight;}
  // destructor
  ~UndirectedEdge() = default;
 private:
  size_t _first, _second;
  double _weight;
};

}  // end namespace algorithms
