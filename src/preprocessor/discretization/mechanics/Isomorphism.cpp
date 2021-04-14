#include "Isomorphism.hpp"
#include <bitset>  // for debugging

namespace discretization {

bool Isomorphism::check(angem::Polyhedron<double> const &p1,
                        angem::Polyhedron<double> const &p2)
{
  // build graph for p1
  Graph const g1 = build_edge_graph_(p1);
  Graph const g2 = build_edge_graph_(p2);
  std::vector<size_t> no_perm(g1.size());
  std::iota(no_perm.begin(), no_perm.end(), 0);
  std::vector<size_t> c1 = compress_(g1, no_perm);
  std::vector<size_t> perm(g2.size());
  std::iota(perm.begin(), perm.end(), 0);

  do {
    auto c2 = compress_(g2, perm);
    if ( c1 == c2 )
    {
      std::cout << "success" << std::endl;
      return true;
    }

  } while (std::next_permutation(perm.begin(), perm.end()));

  std::cout << "fail" << std::endl;
  return false;

  exit(0);
}

std::vector<uint64_t> Isomorphism::compress_(Graph const & g,
                                             std::vector<size_t> const & perm)
{
  /* Must have less vertices than 64! in a polyhedron in order for this to work
     idea: jam colums of a 2D matrix into a single number, e.g.
     Example:
     | 1 0 0 1 0 |
     | 0 0 1 0 1 |  --> [10010, 00101, 11100] (in binary)
     | 1 1 1 0 0 |
  */
  assert( g.size() < 64 );
  std::vector<uint64_t> c(g.size(), 0u);
  for (size_t v = 0; v < g.size(); ++v)
    for (size_t w : g.adj(v))
      c[perm[v]] |= (1 << perm[w]);  // set w'th bit to 1 in row v

  for (auto row : c)
  {
    std::cout << std::bitset<64>(row) << std::endl;
  }

  return c;
}

Graph Isomorphism::build_edge_graph_(angem::Polyhedron<double> const &poly)
{
  algorithms::Graph g(poly.get_points().size());
  for (auto const & edge : poly.get_edges())
  {
    g.add(edge.first, edge.second);
    g.add(edge.second, edge.first);
  }
  return g;
}


}  // end namespace discretization
