#include "Isomorphism.hpp"
#include "angem/utils.hpp"
#include "EdgePath.hpp"
// #include <bitset>  // for debugging

namespace discretization {

using Point = angem::Point<3,double>;
using Polyhedron = angem::Polyhedron<double>;
using std::vector;
using algorithms::Graph;

Isomorphism::Isomorphism(angem::Polyhedron<double> const &p1,
                         angem::Polyhedron<double> const &p2)
    : _p1(p1), _p2(p2)
{
  // numbers of vertices must match
  // if (p1.get_points().size() != p2.get_points().size())
  //   return {false, std::vector<size_t>()};

  // build vertex connectivity graphs
  Graph g1 = build_vertex_graph_(p1);
  // Graph const g2 = build_edge_graph_(p2);
}


Graph Isomorphism::build_vertex_graph_(Polyhedron const &p)
{
  Graph g(p.get_points().size());
  for (auto const & edge : p.get_edges())
    g.add(edge.first, edge.second);

  algorithms::EdgePath path(0, g, p);

  return g;
}

bool Isomorphism::check() const
{
  return false;
}



}  // end namespace discretization
