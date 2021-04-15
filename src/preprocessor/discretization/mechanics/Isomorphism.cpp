#include "Isomorphism.hpp"
#include <bitset>  // for debugging

namespace discretization {

std::pair<bool, std::vector<size_t>>
Isomorphism::check(angem::Polyhedron<double> const &p1,
                   angem::Polyhedron<double> const &p2)
{
  // numbers of vertices must match
  if (p1.get_points().size() != p2.get_points().size())
    return {false, std::vector<size_t>()};

  // build vertex connectivity graphs
  Graph const g1 = build_edge_graph_(p1);
  Graph const g2 = build_edge_graph_(p2);
  // do not renumber vertices in the first graph
  std::vector<size_t> no_perm(g1.size());
  std::iota(no_perm.begin(), no_perm.end(), 0);
  // compress graph for a cheaper comparison
  std::vector<size_t> c1 = compress_(g1, no_perm);

  // keep permuting vertex numbering in the second polyhedron until
  // we get a matching graph or run out of options
  std::vector<size_t> perm(g2.size());
  std::iota(perm.begin(), perm.end(), 0);
  do {
    auto c2 = compress_(g2, perm);
    if ( c1 == c2 )
    {
      std::cout << "success" << std::endl;
      return {true, perm};
    }

  } while (std::next_permutation(perm.begin(), perm.end()));

  std::cout << "fail" << std::endl;
  return {false, std::vector<size_t>()};
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
