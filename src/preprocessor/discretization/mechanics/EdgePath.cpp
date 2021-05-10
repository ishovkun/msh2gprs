#include "EdgePath.hpp"
#include <algorithm>  // all_of
#include <limits>     // numeric_limits

namespace algorithms {

using namespace angem;
using Basis = Basis<3,double>;
using Point = Point<3,double>;
using Polyhedron = Polyhedron<double>;
using Plane = Plane<double>;
static constexpr size_t unmapped = std::numeric_limits<size_t>::max();

EdgePath::EdgePath(size_t source, size_t start_edge,
                   Graph & g, Polyhedron const & poly)
    : EdgePathBase(source, start_edge, g, poly)
{
  size_t const start_edge_global = _g.adj(_s)[_se];
  dfs_(_s, start_edge_global);
  _exists = std::all_of(_visited.begin(), _visited.end(),
      [](bool i) {return i;});
}

// EdgePath::EdgePath(size_t source, size_t start_edge, Graph & g,
//                    angem::Polyhedron<double> const & poly,
//                    EdgePath const & follow)
//     : _g(g), _s(source), _se(start_edge), _verts(poly.get_points()),
//       _path(follow.get_local_edge_path())
// {
//   _basises.resize(_g.nv(), nullptr);
//   _visited.resize(_g.ne(), false);
//   _c = poly.center();
//   // _mapping.resize(_g.nv(), unmapped);
//   size_t const start_edge_global = _g.adj(_s)[_se];
//   _exists = dfs_follow_(_s, start_edge_global, 0);
//   _exists &= std::all_of(_visited.begin(), _visited.end(),
//                          [](bool i) {return i;});
// }

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
      _edge_path.push_back(ie + 1);
      size_t const w = _g.edge(e).other(v);
      dfs_(w, e);
      _vertex_path.push_back(v);
    }

  if (!(incoming == _se && v == _s))
  {
    auto const it = std::find(edge_indices.begin(), edge_indices.end(), incoming);
    size_t const ie = std::distance(edge_indices.begin(), it);
    _edge_path.push_back(-(ie + 1));
  }
}

}  // end namespace algorithms
