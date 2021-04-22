#pragma once
#include "Graph.hpp"
#include "angem/Polyhedron.hpp"
#include "angem/Point.hpp"
#include "angem/Plane.hpp"
#include <memory>

namespace algorithms {

class EdgePath {
  using Basis = angem::Basis<3,double>;
  using Point = angem::Point<3,double>;
  using Polyhedron = angem::Polyhedron<double>;
  using Plane = angem::Plane<double>;

 public:
  // construct path on a polyhedron
  EdgePath(size_t source, size_t start_edge,
           Graph & g, Polyhedron const & poly);

  // see if we can follow the path from a polyhedron
  EdgePath(size_t source, size_t start_edge,
           Graph & g,
           angem::Polyhedron<double> const & poly,
           EdgePath const & follow);

  // returns path of local edge indices (locality : adjacent to the current vertex)
  std::vector<int> const & get_local_edge_path() const {return _path;}
  // returns vertex indices encountered in the order of path (they repeat of course)
  std::vector<size_t> const & get_vertex_path() const {return _vertex_path;}
  bool exist() const {return _exists;}

  virtual ~EdgePath() = default;

 protected:
  void build_basis_(size_t v, size_t edge_to);
  void dfs_(size_t v, size_t edge_to);
  bool dfs_follow_(size_t v, size_t incoming, size_t path_idx);
  void sort_edges_ccw_(size_t v);

  // vertex-vertex graph built on the polyhedron
  Graph _g;
  // coordinates of polyhedron vertices
  std::vector<angem::Point<3,double>> const & _verts;
  // source vertex (first in the path)
  size_t _s{0};
  // first local edge index in the in the path
  size_t _se{0};
  // basis defined in each vertex for identical local edge ordering
  std::vector<std::shared_ptr<Plane>> _basises;
  // center of polyhedron
  angem::Point<3,double> _c;
  // whether an edge has been visited
  std::vector<bool> _visited;
  // path of local edge indices (locality : adjacent to the current vertex)
  std::vector<int> _path;
  // vertex indices encountered in the order of path (they repeat of course)
  std::vector<size_t> _vertex_path;
  // true if able to find path
  bool _exists;
};

}  // end namespace algorithms
