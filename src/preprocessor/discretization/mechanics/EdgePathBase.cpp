#include "EdgePathBase.hpp"

namespace algorithms {

EdgePathBase::EdgePathBase(size_t source, size_t start_edge,
                           Graph & g, Polyhedron const & poly)
    : _s(source), _se(start_edge), _g(g), _verts(poly.get_points()), _center(poly.center())
{
  _basises.resize(_g.nv(), nullptr);
  _visited.resize(_g.ne(), false);
}

void EdgePathBase::build_basis_(size_t v, size_t edge_to)
{
  auto normal = (_verts[v] - _center).normalize();
  _basises[v] = std::make_shared<Plane>(_center, normal);
  size_t const w = _g.edge(edge_to).other(v);
  auto p = _basises[v]->project_point(_verts[w]);
  auto t1 = (p - _center).normalize();
  auto t2 = normal.cross(t1).normalize();
  _basises[v]->set_basis(Basis(t1, t2, normal));
}

void EdgePathBase::sort_edges_ccw_(size_t v)
{
  std::vector<Point> neighbor_vertices(_g.degree(v));
  auto const plane = _basises[v];
  auto edges = _g.adj(v);
  for (size_t i = 0; i < _g.degree(v); ++i)
    neighbor_vertices[i] = _verts[ _g.edge(edges[i]).other(v) ];
  auto order = angem::Polygon<double>::order_ccw(neighbor_vertices, *plane);
  angem::reorder_from(_g.adj(v), order);
}

}  // end namespace algorithms
