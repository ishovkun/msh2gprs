#pragma once
#include "angem/Polyhedron.hpp"
#include "Graph.hpp"

namespace discretization
{

/*
 * Check if two polyhedra are combinatorially equivalent (isomorphic).
 * Two polyhedra are called combinatorially equivalent if their vertex-vertex graphs
 * are isomorphic.
 * Yost, D. (2007). Some indecomposable polyhedra. Optimization, 56(5-6):715–724.
 *
 * This class implements an isomorphism algorithm insipired by the paper
 * Sugihara, K. (1984). An n log n algorithm for determining the congruity of polyhedra.
 * Journal of Computer and System Sciences, 29(1):36–47.
 *
 * The idea is to distinguish edge paths in a polyhedron geometrically (e.g. turn left/right).
 * We first build a primary path on the reference polyhedron, and try to repeat it
 * over the edges of the second polyhedron. If we can do it and visit all the edges,
 * then two polyhedra are isomorphic.
 *
 * Additionally, upon a positive test, this class can return the mapping of the vertices of
 * the current polyhedron to those in the reference polyhedron.
 *
 * The time complexity of the check is O(n_vertices^2).
 * Reduces to O(n) if the initial vertex ordering is correct.
 */
class Isomorphism {
 public:
  /* Constructor. Takes references for two polyhedra to be tested. */
  Isomorphism(angem::Polyhedron<double> const &reference,
              angem::Polyhedron<double> const &current);
  // returns true if polyhedras are isomorphic
  bool check() const {return !_ordering.empty();}
  // returns the mapping of the vertices of the current polyhedron to those in the reference polyhedron.
  std::vector<size_t> const & ordering() const {return _ordering;}
  /* Destructor */
  virtual ~Isomorphism() = default;

 private:
  algorithms::Graph build_vertex_graph_(angem::Polyhedron<double> const &poly);

  angem::Polyhedron<double> const & _p1;
  angem::Polyhedron<double> const & _p2;
  // vertex ordering built upon concluding equivalence
  std::vector<size_t> _ordering;
};

}
