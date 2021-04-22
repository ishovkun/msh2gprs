#include "EdgePath.hpp"
#include <algorithm>  // all_of

namespace algorithms {

using namespace angem;
using Basis = Basis<3,double>;
using Point = Point<3,double>;
using Polyhedron = Polyhedron<double>;
using Plane = Plane<double>;

EdgePath::EdgePath(size_t source, size_t start_edge,
                   Graph & g, Polyhedron const & poly)
    : _s(source), _se(start_edge), _g(g), _verts(poly.get_points())
{
  _basises.resize(_g.nv(), nullptr);
  _visited.resize(_g.ne(), false);
  _c = poly.center();

  // initialize basis in source vertex starting from edge se
  auto const * edge = g.adj(_s)[_se];
  build_basis_(_s, *edge);
  sort_edges_ccw_(_s);
  dfs_(_s);

  _exists = std::all_of(_visited.begin(), _visited.end(),
      [](bool i) {return i;});
}

EdgePath::EdgePath(size_t source, size_t start_edge, Graph & g,
                   angem::Polyhedron<double> const & poly,
                   EdgePath const & follow)
    : _g(g), _se(start_edge), _verts(poly.get_points()), _path(follow.get())
{
  _basises.resize(_g.nv(), nullptr);
  _visited.resize(_g.ne(), false);
  _c = poly.center();

  // initialize basis in source vertex starting from edge se
  auto const * edge = g.adj(_s)[_se];
  build_basis_(_s, *edge);
  sort_edges_ccw_(_s);
  dfs_follow_(_s, 0);

  _exists = std::all_of(_visited.begin(), _visited.end(),
      [](bool i) {return i;});
}

void EdgePath::dfs_(size_t v)
{
  std::cout << "visiting vertex " << v << " ";
  auto edges = _g.adj_idx(v);
  for (size_t ie = 0; ie < edges.size(); ++ie) {
    size_t const e = edges[ie];

    if (!_visited[e]) {
      _visited[e] = true;
      _path.push_back(ie+1);
      auto const & edge = _g.edge(e);
      size_t const w = edge.other(v);
      if (!_basises[w]) {
        build_basis_(w, edge);
        sort_edges_ccw_(w);
      }
      std::cout << "following edge " << ie + 1 << " (" << e << ") "
                << "to vertex " << w;
      std::cout << std::endl;
      dfs_(w);
      auto new_adj = _g.adj_idx(w);
      auto loc = std::find(new_adj.begin(), new_adj.end(), e);
      int back_edge = std::distance(new_adj.begin(), loc);
      _path.push_back(-(back_edge + 1));
    }
  }
}

bool EdgePath::dfs_follow_(size_t v, size_t cur)
{
  std::cout << "visiting vertex " << v << " ";
  if (cur >= _path.size()) return true;

  int const ie = (_path[cur] > 0) ? _path[cur] - 1 : _path[cur] + 1;

  if (std::abs(ie) > _g.degree(v)) {
    std::cout << "edge not present " << ie << std::endl;
    return false;
  }

  auto edges = _g.adj_idx(v);
  size_t e = edges[std::abs(ie)];
  auto const & edge = _g.edge(e);

  std::cout << "following edge " << _path[cur] << " (" << e << ") "
            << "to vertex " << edge.other(v) << " ";

  if (ie >= 0 && _visited[e])
  {
    std::cout << "already visited" << std::endl;
    return false;
  }
  else if (ie < 0 && !_visited[e])
  {
    std::cout << "backward edge not visited" << std::endl;
    return false;
  }
  std::cout << std::endl;

  _visited[e] = true;

  size_t const w = edge.other(v);
  if (!_basises[w]) {
    build_basis_(w, edge);
    sort_edges_ccw_(w);
  }

  return dfs_follow_(w, cur+1);
}

void EdgePath::build_basis_(size_t v, Edge const & incoming)
{
  auto normal = (_verts[v] - _c).normalize();
  _basises[v] = std::make_shared<Plane>(_c, normal);
  size_t w = incoming.other(v);
  auto p = _basises[v]->project_point(_verts[w]);
  auto t1 = (p - _c).normalize();
  auto t2 = normal.cross(t1).normalize();
  _basises[v]->set_basis(Basis(t1, t2, normal));
}

void EdgePath::sort_edges_ccw_(size_t v)
{
  std::vector<Point> neighbor_vertices(_g.degree(v));
  auto const plane = _basises[v];
  auto edges = _g.adj(v);
  for (size_t i = 0; i < _g.degree(v); ++i)
    neighbor_vertices[i] = _verts[ edges[i]->other(v) ];
  auto order = angem::Polygon<double>::order_ccw(neighbor_vertices, *plane);
  // std::cout << "order ";
  // for (auto i : order) std::cout << i << " ";
  // std::cout << std::endl;
  // std::cout << "\npre-order ";
  // for (auto i : _g.adj_idx(v)) std::cout << i << " ";
  // std::cout << std::endl;
  _g.reorder_vertex_edges(v, order);
  // std::cout << "post-order ";
  // for (auto i : _g.adj_idx(v)) std::cout << i << " ";
  // std::cout << std::endl;
}


}  // end namespace algorithms
