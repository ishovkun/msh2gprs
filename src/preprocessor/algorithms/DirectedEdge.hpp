#pragma once
#include <cstddef>  // size_t
#include <limits>   // numeric_limits

namespace algorithms {

class DirectedEdge {
 private:
  static constexpr size_t undefined = std::numeric_limits<size_t>::max();
  static constexpr double inf = std::numeric_limits<double>::max();

 public:
  DirectedEdge(const size_t v, const size_t w, const double weight) noexcept
      : _v(v), _w(w), _weight(weight)
  {}

  DirectedEdge() noexcept
      : _v(undefined), _w(undefined), _weight(inf)
  {}

  DirectedEdge & operator=(const DirectedEdge & other)
  {
    _v = other.from();
    _w = other.to();
    _weight = other.weight();
    return *this;
  }

  inline const size_t from() const noexcept {return _v;}
  inline const size_t to() const noexcept {return _w;}
  inline const double weight() const noexcept {return _weight;}
  inline void set_weight(double new_weight) noexcept {_weight = new_weight;}

  virtual ~DirectedEdge() = default;

 private:
  size_t _v, _w;
  double _weight;
};


}  // end namespace algorithms
