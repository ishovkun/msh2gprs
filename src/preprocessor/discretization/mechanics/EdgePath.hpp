#pragma once
#include "Graph.hpp"
#include "angem/Polyhedron.hpp"
#include "angem/Point.hpp"
#include "angem/Basis.hpp"
#include <memory>

namespace algorithms {

class EdgePath {
  using Basis = angem::Basis<3,double>;
  using Point = angem::Point<3,double>;
  using Polyhedron = angem::Polyhedron<double>;

 public:
  EdgePath(size_t source, Graph & g, Polyhedron const & poly);

  std::vector<size_t> & get_();

  virtual ~EdgePath() = default;

 private:
  void build_basis_(size_t v, Point dir_in);
  void dfs_(size_t v);

  Graph & _g;
  std::vector<angem::Point<3,double>> const & _verts;  // coord
  size_t _s;  // source
  std::vector<std::shared_ptr<Basis>> _basises;  // should be unique but I'm pissed
  angem::Point<3,double> _c;
  std::vector<bool> _visited;
  std::vector<size_t> _path;
};

}  // end namespace algorithms
