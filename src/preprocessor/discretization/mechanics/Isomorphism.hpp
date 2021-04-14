#pragma once
#include "angem/Polyhedron.hpp"
#include "Graph.hpp"

namespace discretization
{

using algorithms::Graph;

class Isomorphism {
 public:
  Isomorphism() = delete;
  virtual ~Isomorphism() = delete;

  static bool check(angem::Polyhedron<double> const &p1,
                    angem::Polyhedron<double> const &p2);

 protected:
  static Graph build_edge_graph_(angem::Polyhedron<double> const &poly);
  static std::vector<uint64_t> compress_(Graph const & g,
                                         std::vector<size_t> const & perm);
};

}
