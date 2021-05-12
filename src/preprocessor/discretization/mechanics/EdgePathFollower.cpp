#include "EdgePathFollower.hpp"

static constexpr size_t unmapped = std::numeric_limits<size_t>::max();

namespace algorithms {

EdgePathFollower::EdgePathFollower(size_t source, size_t start_edge, Graph & g,
                                   angem::Polyhedron<double> const & poly,
                                   EdgePath const & follow)
    : EdgePathBase(source, start_edge, g, poly),
      _edge_path(follow.get_local_edge_path()),
      _base_vertex_path(follow.get_vertex_path())
{
   // _vertex_path.assign(_base_vertex_path.size(), unmapped);
   _mapping.assign(g.nv(), unmapped);

  size_t const start_edge_global = _g.adj(_s)[_se];
  _exists = dfs_(_s, start_edge_global, 0);
  // for (auto v : _vertex_path)
  //   std::cout << v << " ";
  // std::cout << std::endl;
  // std::cout << std::endl;

  // make sure we visited all the edges
  _exists &= std::all_of(_visited.begin(), _visited.end(), [](bool i) {return i;});
}

bool EdgePathFollower::dfs_(size_t v, size_t incoming, size_t path_idx)
{
  if (_mapping[v] != unmapped &&  _mapping[v] != _base_vertex_path[path_idx])
  {
    // std::cout << "mapping violation:" << std::endl;
    // std::cout << "visiting vertex " << v;
    // std::cout << " mapped to " << _mapping[v] << std::endl;
    // std::cout << "but current match says " << _base_vertex_path[path_idx] << std::endl;
    // std::cout << "path idx = " << path_idx << " (path size = " << _base_vertex_path.size() << ")" << std::endl;

    return false;
  }

  _mapping[v] = _base_vertex_path[path_idx];
  _vertex_path.push_back(v);

  if (path_idx == _edge_path.size() - 1)  // base case
  {
    // std::cout << "finish: OK" << std::endl;
    return true;
  }

  if (!_basises[v])
  {
    build_basis_(v, incoming);
    sort_edges_ccw_(v);
  }

  auto const & edge_indices = _g.adj(v);
  int const step = _edge_path[path_idx];  // negative - we are going back :-)
  bool const goback = step < 0;
  int const ie = abs(step) - 1;  // local_edge_index

  if (ie >= _g.degree(v))  // cannot proceed along this edge since this vertex has less edges
  {
    // std::cout << "degree not match" << std::endl;
    return false;
  }

  size_t const e = edge_indices[ie];
  // this edges already visited and we are not returning, should not happen for primary path
  if (_visited[e] && !goback)
  {
    // std::cout << "edge visited and not going back" << std::endl;
    return false;
  }

  _visited[e] = true;
  size_t const w = _g.edge(e).other(v);  // next vertex along the edge e
  return dfs_(w, e, path_idx + 1);
}

}  // end namespace algorithms
