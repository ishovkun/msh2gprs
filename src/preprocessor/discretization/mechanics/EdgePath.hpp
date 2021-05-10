#pragma once
#include "EdgePathBase.hpp"
// #include "Graph.hpp"
// #include "angem/Polyhedron.hpp"
// #include "angem/Point.hpp"
// #include "angem/Plane.hpp"
#include <memory>

namespace algorithms {

class EdgePath : public EdgePathBase {

  using Basis = angem::Basis<3,double>;
  using Point = angem::Point<3,double>;
  using Polyhedron = angem::Polyhedron<double>;
  using Plane = angem::Plane<double>;

 public:
  // construct path on a polyhedron
  EdgePath(size_t source, size_t start_edge,
           Graph & g, Polyhedron const & poly);

  // see if we can follow the path from a polyhedron
  // EdgePath(size_t source, size_t start_edge,
  //          Graph & g,
  //          angem::Polyhedron<double> const & poly,
  //          EdgePath const & follow);

  // returns path of local edge indices (locality : adjacent to the current vertex)
  std::vector<int> const & get_local_edge_path() const {return _edge_path;}
  // returns vertex indices encountered in the order of path (they repeat of course)
  std::vector<size_t> const & get_vertex_path() const {return _vertex_path;}
  // returns true if was able to construct primary path
  bool exist() const {return _exists;}

  virtual ~EdgePath() = default;

 protected:
  void dfs_(size_t v, size_t edge_to);

  // path of local edge indices (locality : adjacent to the current vertex)
  std::vector<int> _edge_path;
  // vertex indices encountered in the order of path (they repeat of course)
  std::vector<size_t> _vertex_path;
  // true if able to find path
  bool _exists;
  // // vertex mapping for the following path
  // std::vector<size_t> _mapping;
};

}  // end namespace algorithms
