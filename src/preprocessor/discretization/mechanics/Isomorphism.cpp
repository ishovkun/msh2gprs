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
  // build vertex connectivity graphs
  Graph g1 = build_vertex_graph_(p1);
  // std::cout << "\nbuilding path 1" << std::endl;
  algorithms::EdgePath path1(/*source */ 0, /*edge */ 0, g1, p1);
  // std::cout << "\nbuilding path 3" << std::endl;
  // algorithms::EdgePath path3(/*source */ 0, /*edge */ 1, g1, p1);
  // assert( path1.exists() );
  // exit(0);
  //
  std::cout << "found path: ";
  for ( auto e : path1.get() )
    std::cout << e << " ";
  std::cout << "(total edges = " << g1.ne() << ")"<< std::endl;

  // std::cout << "testing: " << std::endl;
  // for (size_t ie = 0; ie < g1.degree(0); ++ie) {
  //   // std::cout << "\nstarting edge = " << ie << std::endl;
  //   algorithms::EdgePath path2(0, 0, g1, p1, path1);
  //   // std::cout << "exists? " << path2.exist() << std::endl;
  // }

  // exit(0);
  //
  Graph g2 = build_vertex_graph_(p2);
  for (size_t start = 0; start < g2.nv(); ++start)
    for (size_t ie = 0; ie < g2.degree(start); ++ie)
    {
      std::cout << "\n=============================\n"
                << "trying numeration v = " << start << " e = " << ie << std::endl;
      algorithms::EdgePath path2(start, ie, g2, p2, path1);
      // algorithms::EdgePath path2(start, ie, g1, p1, path1);
      if (path2.exist())
      {
        std::cout << "fuck yeah" << std::endl;
        exit(0);
        break;
      }
      // else std::cout << "fuck no" << std::endl;
      // exit(0);
    }
  std::cout << "fuck no" << std::endl;

  exit(0);
}


Graph Isomorphism::build_vertex_graph_(Polyhedron const &p)
{
  Graph g(p.get_points().size());
  for (auto const & edge : p.get_edges())
    g.add(edge.first, edge.second);

  return g;
}

bool Isomorphism::check() const
{
  return false;
}



}  // end namespace discretization
