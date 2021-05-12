#include "Isomorphism.hpp"
#include "angem/utils.hpp"
#include "EdgePath.hpp"
#include "EdgePathFollower.hpp"
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
  // build vertex connectivity graphs
  Graph g1 = build_vertex_graph_(p1);
  algorithms::EdgePath path1(/*source */ 0, /*edge */ 0, g1, p1);
  if ( !path1.exist() ) throw std::runtime_error("Cannot build primary path on the reference element");

  // std::cout << "vertex path 1: ";
  // for (auto v : path1.get_vertex_path())
  //   std::cout << v << " ";
  // std::cout << " ("<< path1.get_vertex_path().size() << ")" << std::endl;
  // std::cout << "edge path 1 : ";
  // for (auto e : path1.get_local_edge_path())
  //   std::cout << e << " ";
  // std::cout << " (" << path1.get_local_edge_path().size() << ")" << std::endl;

  Graph g2 = build_vertex_graph_(p2);
  for (size_t start = 0; start < g2.nv(); ++start) {
    for (size_t ie = 0; ie < g2.degree(start); ++ie) {
      algorithms::EdgePathFollower path2(start, ie, g2, p2, path1);

      // std::cout << "vp2: ";
      // for (auto v : path2.get_vertex_path())
      //   std::cout << v << " ";
      // std::cout << " ("<< path2.get_vertex_path().size() << ")" << std::endl;
      // for (auto e : path2.get_local_edge_path())
      //   std::cout << e << " ";
      // std::cout << " (" << path2.get_local_edge_path().size() << ")" << std::endl;

      if (path2.exist()) {
        _ordering = path2.get_vertex_mapping();
        return;
      }
    }
  }
}


Graph Isomorphism::build_vertex_graph_(Polyhedron const &p)
{
  Graph g(p.get_points().size());
  for (auto const & edge : p.get_edges())
    g.add(edge.first, edge.second);
  return g;
}

}  // end namespace discretization
