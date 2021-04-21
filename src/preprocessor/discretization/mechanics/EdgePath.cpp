#include "EdgePath.hpp"

namespace algorithms {

using Basis = angem::Basis<3,double>;
using Point = angem::Point<3,double>;
using Polyhedron = angem::Polyhedron<double>;

EdgePath::EdgePath(size_t source, Graph & g, Polyhedron const & poly)
    : _s(source), _g(g), _verts(poly.get_points())
{
  _basises.resize(_g.nv(), nullptr);
  _visited.resize(_g.ne(), false);

  _c = poly.center();
  auto const * edge = *_g.adj(_s).begin();
  auto dir1 = (_verts[edge->other(_s)] - _verts[_s]).normalize();
  build_basis_(_s, dir1);
  dfs_(_s);
}

void EdgePath::dfs_(size_t v)
{
  for (size_t e : _g.adj_idx(v))
    if (!_visited[e])
    {
      _visited[e] = true;
      auto const & edge = _g.edge(e);
      size_t const w = edge.other(v);
      if (!_basises[w])
        build_basis_(v, _verts[w] - _verts[v]);
      // TODO: sort outgoing edges ccw wrt basis
      throw std::runtime_error("not implemented");
    }
  // _visited[v] = true;
  // for (const size_t w: _gr.adj(v))
  //   if (!_marked[w])
  //   {
  //     _edge[w] = v;
  //     dfs(w);
  //   }
}

void EdgePath::build_basis_(size_t v, EdgePath::Point dir1)
{
  auto dir2 = (_c - _verts[_s]).normalize();
  auto dir3 = dir1.cross(dir2);
  _basises[v] = std::make_shared<Basis>(dir1, dir2, dir3);
}

}  // end namespace algorithms
