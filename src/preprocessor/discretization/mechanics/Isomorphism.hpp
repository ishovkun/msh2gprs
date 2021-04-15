#pragma once
#include "angem/Polyhedron.hpp"
#include "Graph.hpp"

namespace discretization
{

using algorithms::Graph;

class Isomorphism {
 public:
  // check if two polyhedral are combinatorically equivalent (isomorphic)
  // and, if so, return vertex renumbering for the second polyhedron
  static std::pair<bool, std::vector<size_t>>
  check(angem::Polyhedron<double> const &master,
        angem::Polyhedron<double> const &current);

 private:
  Isomorphism() = delete;
  virtual ~Isomorphism() = delete;
  static Graph build_edge_graph_(angem::Polyhedron<double> const &poly);
  static std::vector<uint64_t> compress_(Graph const & g,
                                         std::vector<size_t> const & perm);
};

}
