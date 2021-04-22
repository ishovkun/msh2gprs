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

  size_t start_edge_global = _g.adj(_s)[_se];
  dfs_(_s, start_edge_global);

  _exists = std::all_of(_visited.begin(), _visited.end(),
      [](bool i) {return i;});
}

EdgePath::EdgePath(size_t source, size_t start_edge, Graph & g,
                   angem::Polyhedron<double> const & poly,
                   EdgePath const & follow)
    : _g(g), _s(source), _se(start_edge), _verts(poly.get_points()),
      _path(follow.get_local_edge_path())
{
  _basises.resize(_g.nv(), nullptr);
  _visited.resize(_g.ne(), false);
  _c = poly.center();
  size_t const start_edge_global = _g.adj(_s)[_se];
  _exists = dfs_follow_(_s, start_edge_global, 0);
  _exists &= std::all_of(_visited.begin(), _visited.end(),
                         [](bool i) {return i;});
}

// for debug
void print_angles(std::vector<Point> const & points,
                  Plane const & plane, double eps)
{
  size_t const np = points.size();
  std::vector<Point> local(np);
  for (size_t i = 0; i < np; ++i)
    local[i] = plane.local_coordinates(points[i]);

  std::vector<double> angles(np, 0);
  for (size_t i = 0; i < np; ++i) {
    angles[i] = std::atan2(local[i][1], local[i][0]);
    if (angles[i] < 0 && angles[i] > -eps)
      angles[i] = 0;
    if (angles[i] < 0)
      angles[i] = 2*M_PI - std::fabs(angles[i]);
    std::cout << "angles[" << i << "] = " << angles[i] << std::endl;
  }

}

void EdgePath::dfs_(size_t v, size_t incoming)
{
  // std::cout << "visiting ";
  // if (!(!_visited[incoming] && v == _s))  // skip first
  //   std::cout << _g.edge(incoming).other(v) << "-" << v;
  // else std::cout << v;
  // std::cout << std::endl;

  _vertex_path.push_back(v);

  if (!_basises[v])
  {
    build_basis_(v, incoming);
    sort_edges_ccw_(v);
  }

  auto const edge_indices = _g.adj(v);
  for (size_t ie = 0; ie < _g.degree(v); ++ie)
    if (!_visited[edge_indices[ie]])
    {
      size_t const e = edge_indices[ie];
      _visited[e] = true;
      _path.push_back(ie + 1);
      size_t const w = _g.edge(e).other(v);
      dfs_(w, e);
    }

  if (!(incoming == _se && v == _s))
  {
    auto const it = std::find(edge_indices.begin(), edge_indices.end(), incoming);
    size_t const ie = std::distance(edge_indices.begin(), it);
    _path.push_back(-(ie + 1));
  }
}

bool EdgePath::dfs_follow_(size_t v, size_t incoming, size_t path_idx)
{
  if (path_idx == 0 || (path_idx > 0 && _path[path_idx-1] >= 0))
    _vertex_path.push_back(v);

  if (path_idx == _path.size())  // base case
    return true;

  // if (!(!_visited[incoming] && v == _s))  // skip first
  //   std::cout << "visiting " << _g.edge(incoming).other(v) << "-" << v << std::endl;
  // else std::cout << "visiting " << v << std::endl;

  if (!_basises[v])
  {
    build_basis_(v, incoming);
    sort_edges_ccw_(v);
  }

  auto const & edge_indices = _g.adj(v);
  int step = _path[path_idx];
  bool goback = step < 0;
  int ie = abs(step) - 1;  // local_edge_index
  // std::cout << "follow " << step << " (" << edge_indices[ie] << ") "
  //           << "(path " <<  path_idx
  //           << " out of " << _path.size() << ")"
  //           << std::endl;

  if (ie >= _g.degree(v)) {
    // std::cout << "invalid edge "  << std::endl;
    return false;
  }

  size_t e = edge_indices[ie];
  if (_visited[e] && !goback) {
    // std::cout << "already visited" << std::endl;
    return false;
  }

  // else if (_visited[e] && goback)
  _visited[e] = true;
  size_t w = _g.edge(e).other(v);
  return dfs_follow_(w, e, path_idx + 1);
}

void EdgePath::build_basis_(size_t v, size_t edge_to)
{
  auto normal = (_verts[v] - _c).normalize();
  _basises[v] = std::make_shared<Plane>(_c, normal);
  size_t const w = _g.edge(edge_to).other(v);
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
    neighbor_vertices[i] = _verts[ _g.edge(edges[i]).other(v) ];
  auto order = angem::Polygon<double>::order_ccw(neighbor_vertices, *plane);
  // std::cout << "order ";
  // for (auto i : order) std::cout << i << " ";
  // std::cout << std::endl;
  // std::cout << "\npre-order ";
  // for (auto i : _g.adj_idx(v)) std::cout << i << " ";
  // std::cout << std::endl;
  // std::cout << "order ";
  // for (auto i : order) std::cout << i <<  " ";
  // std::cout << std::endl;
  // std::vector<Point> pts;
  // for (auto e : _g.adj_idx(v))
  //   pts.push_back( _verts[_g.edge(e).other(v)] );
  // print_angles(pts, *_basises[v], 1e-8);

  angem::reorder_from(_g.adj(v), order);

  // std::cout << "post-order ";
  // for (auto i : _g.adj_idx(v)) std::cout << i << " ";
  // std::cout << std::endl;
  // pts.clear();
  // for (auto e : _g.adj_idx(v))
  //   pts.push_back( _verts[_g.edge(e).other(v)] );
  // print_angles(pts, *_basises[v], 1e-8);

}


}  // end namespace algorithms
