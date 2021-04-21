#include "Isomorphism.hpp"
#include "angem/utils.hpp"
// #include <bitset>  // for debugging

namespace discretization {

using Point = angem::Point<3,double>;
using Polyhedron = angem::Polyhedron<double>;
using std::vector;

inline double dot_product_plane_basis_normal(angem::Point<3,double> const & t1,
                                             angem::Point<3,double> const & t2,
                                             angem::Point<3,double> const & n)
{
  return t1.cross(t2).dot(n) > 0;
}

vector<size_t> uncompress(std::vector<uint64_t> const & c, size_t row)
{
  size_t const nv = c.size();
  vector<size_t> uc;
  size_t x = c[row];
  for (size_t i = 0; i < nv; ++i)
    if ( x & (1 << i) ) uc.push_back(i);
  return uc;
}

bool check_orientation(Polyhedron const & p1,
                       Polyhedron const & p2,
                       std::vector<size_t> order,
                       std::vector<uint64_t> const & c)
{
  auto const & v1 = p1.get_points();
  auto v2 = p2.get_points();
  angem::reorder(v2, order);
  std::vector<size_t> adj = uncompress(c, 0);
  double const d1 = dot_product_plane_basis_normal(v1[adj[0]] - v1[0],
                                                   v1[adj[1]] - v1[0],
                                                   v1[adj[2]] - v1[0]);
  double const d2 = dot_product_plane_basis_normal(v2[adj[0]] - v2[0],
                                                   v2[adj[1]] - v2[0],
                                                   v2[adj[2]] - v2[0]);
  return d1 * d2 > 0;
}

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
    if ( c1 == c2 && check_orientation(p1, p2, perm, c1))
    {
      return {true, perm};
    }

  } while (std::next_permutation(perm.begin(), perm.end()));

  return {false, std::vector<size_t>()};
}

std::vector<uint64_t> Isomorphism::compress_(Graph const & g,
                                             std::vector<size_t> const & perm)
{
  /*Idea: jam colums of a 2D matrix into a single number, e.g.
    NOTE: Must have less vertices than 64 in a polyhedron in order for this to work
    Example:
    | 1 0 0 1 0 |
    | 0 0 1 0 1 |  --> [10010, 00101, 11100] (in binary)
    | 1 1 1 0 0 |
  */
  assert( g.size() < 64 && "Cannot compress polyhedrons with >= 64 vertices");

  std::vector<uint64_t> c(g.size(), 0u);
  for (size_t v = 0; v < g.size(); ++v)
    for (size_t w : g.adj(v))
      c[perm[v]] |= (1 << perm[w]);  // set w'th bit to 1 in row v

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
