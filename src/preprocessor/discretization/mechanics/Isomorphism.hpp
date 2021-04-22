#pragma once
#include "angem/Polyhedron.hpp"
#include "Graph.hpp"

namespace discretization
{

class Isomorphism {
 public:
  Isomorphism(angem::Polyhedron<double> const &master,
              angem::Polyhedron<double> const &current);
  // returns true if polyhedras are topologically equivalent
  bool check() const {return !_ordering.empty();}
  std::vector<size_t> const & ordering() const {return _ordering;}
  virtual ~Isomorphism() = default;

 private:
  algorithms::Graph build_vertex_graph_(angem::Polyhedron<double> const &poly);
  void build_ordering_(std::vector<size_t> const & path1,
                       std::vector<size_t> const & path2);

  angem::Polyhedron<double> const & _p1;
  angem::Polyhedron<double> const & _p2;
  // vertex ordering built upon concluding equivalence
  std::vector<size_t> _ordering;
};

}
