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

  std::vector<int> const & get() const {return _path;}
  bool exist() const {return _exists;}

  virtual ~EdgePath() = default;

 protected:
  void build_basis_(size_t v, Edge const & first_edge);
  void dfs_(size_t v);
  bool dfs_follow_(size_t v, size_t cur_edge);
  void sort_edges_ccw_(size_t v);

  Graph _g;
  std::vector<angem::Point<3,double>> const & _verts;  // coord
  size_t _s{0};  // source
  size_t _se{0};
  std::vector<std::shared_ptr<Plane>> _basises;  // should be unique but I'm pissed
  angem::Point<3,double> _c;
  std::vector<bool> _visited;
  std::vector<int> _path;
  bool _exists;
};

}  // end namespace algorithms
