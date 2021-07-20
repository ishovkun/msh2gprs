#pragma once
#include "Graph.hpp"
#include "angem/Polyhedron.hpp"
#include "angem/Point.hpp"
#include "angem/Plane.hpp"
#include <memory>

namespace algorithms {

/* This is a base dummy class to avoid code duplication */
class EdgePathBase {
  using Basis = angem::Basis<3,double>;
  using Point = angem::Point<3,double>;
  using Polyhedron = angem::Polyhedron<double>;
  using Plane = angem::Plane<double>;

 public:
  EdgePathBase(size_t source, size_t start_edge,
               Graph & g, Polyhedron const & poly);

  virtual ~EdgePathBase() = default;

 protected:
  void build_basis_(size_t v, size_t edge_to);
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
  angem::Point<3,double> _center;
  // whether an edge has been visited
  std::vector<bool> _visited;
};

}  // end namespace algorithms
