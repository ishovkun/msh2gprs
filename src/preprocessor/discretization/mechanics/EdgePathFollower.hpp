#pragma once
#include "EdgePath.hpp"

namespace algorithms {

/* This class tries to follow the path built on another polyhedra.
 * If it is successful, all edges are visited, and there is no conflict in
 * vertex-vertex mapping, then we conclude isomorphism */
class EdgePathFollower : public EdgePathBase {
 public:
  EdgePathFollower(size_t source, size_t start_edge, Graph & g,
                   angem::Polyhedron<double> const & poly,
                   EdgePath const & follow);
  // destructor
  virtual ~EdgePathFollower() = default;
  // returns path of local edge indices (locality : adjacent to the current vertex)
  std::vector<int> const & get_local_edge_path() const {return _edge_path;}
  // returns vertex indices encountered in the order of path (they repeat of course)
  std::vector<size_t> const & get_vertex_path() const {return _vertex_path;}
  // returns true if was able to follow the given path
  bool exist() const {return _exists;}
  // return mapping of current path vertices to the base path vertices
  std::vector<size_t> const & get_vertex_mapping() const {return _mapping;}

 protected:
  bool dfs_(size_t v, size_t incoming, size_t path_idx);

  // path of local edge indices (locality : adjacent to the current vertex)
  std::vector<int> _edge_path;
  // vertex indices encountered in the order of path (they repeat of course)
  std::vector<size_t> _base_vertex_path;
  // mapping of the element vertices to the base path
  std::vector<size_t> _mapping;
  // vertex indices encountered in the order of path (they repeat of course)
  std::vector<size_t> _vertex_path;
  // true if able to find path
  bool _exists;
};

}  // end namespace algorithms
