#pragma once
#include <cstddef>  // size_t
#include <cassert>
#include <ostream>
#include <stdexcept>

namespace algorithms {

class FlowEdge {
 public:
  FlowEdge(const size_t v, const size_t w, const double capacity)
      : _from(v), _to(w), _capacity(capacity), _flow(0)
  {}

  virtual ~FlowEdge() = default;
  inline size_t from() const noexcept {return _from;}
  inline size_t to() const noexcept {return _to;}
  inline size_t other(const size_t v) const {
    assert( from() == v || to() == v );
    return (from() == v) ? to() : from() ;
  }

  inline double capacity() const noexcept {return _capacity;}
  inline double flow() const noexcept {return _flow;}
  double residualCapacityTo(const size_t v) const
  {
    if   (v == from()) return flow();
    else if (v == to()) return capacity() - flow();
    else throw std::invalid_argument("vertex" + std::to_string(v) + " is not in edge");
  }

  void addResidualFlowTo(const size_t v, const double delta)
  {
    if   (v == from()) _flow -= delta;
    else if (v == to()) _flow += delta;
    else throw std::invalid_argument("vertex" + std::to_string(v) + " is not in edge");
  }

 private:
  size_t _from, _to;
  double _capacity, _flow;
};

}  // end namespace algorithms
