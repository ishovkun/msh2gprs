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
  bool check() const;
  // TODO: writeme
  std::vector<size_t> const & ordering() const;
  virtual ~Isomorphism() = default;

 private:
  algorithms::Graph build_vertex_graph_(angem::Polyhedron<double> const &poly);

  // static std::vector<uint64_t> compress_(Graph const & g,
  //                                        std::vector<size_t> const & perm);
  angem::Polyhedron<double> const & _p1;
  angem::Polyhedron<double> const & _p2;
};

}
