#ifdef WITH_EIGEN
#include "IntegrationRuleVerticesAverage.hpp"
#include "../FeValues.hpp"

namespace discretization {

using Point = angem::Point<3,double>;

IntegrationRuleVerticesAverage::IntegrationRuleVerticesAverage(PolyhedralElementBase & element)
    : _element(element)
{
  if (_element._face_domains.empty())
    _element._face_domains = _element.create_face_domains_();

  build_tributary_shapes_cells_();
  const auto face_polygons = _element._parent_cell.polyhedron()->get_face_polygons();
  const size_t n_faces = face_polygons.size();
  for (size_t iface=0; iface < n_faces; ++iface)
    build_tributary_shapes_face_(iface, face_polygons[iface]);

  setup_storage_();
  compute_cell_fe_quantities_();
  for (size_t iface=0; iface < n_faces; ++iface)
    compute_face_fe_quantities_(iface);
}

void IntegrationRuleVerticesAverage::setup_storage_()
{
  auto & cell_data = _element._cell_data;
  const auto & basis_functions = _element._basis_functions;
  cell_data.points.resize( _tributary3d.size() );
  for (auto & point : cell_data.points)
  {
    point.values.resize( basis_functions.size(), 0 );
    point.grads.resize( basis_functions.size() );
  }

  cell_data.center.values.resize( basis_functions.size(), 0 );
  cell_data.center.grads.resize( basis_functions.size() );

  const std::vector<const mesh::Face*> parent_faces = _element._parent_cell.faces();
  _element._face_data.resize( _tributary2d.size() );
  for (size_t parent_face = 0; parent_face < _tributary2d.size(); ++parent_face)
  {
    auto & data = _element._face_data[parent_face];
    const size_t n_vertices = parent_faces[parent_face]->vertices().size();
    data.points.resize( _tributary2d[parent_face].size() );
    for (size_t q=0; q<data.points.size(); ++q)
    {
      data.points[q].values.resize(n_vertices);
      data.points[q].grads.resize(n_vertices);
    }
    data.center.values.resize(n_vertices);
    data.center.grads.resize(n_vertices);
  }
}

std::vector<mesh::Edge> get_edges_with_vertex(const size_t vertex,
                                              const mesh::Face & face)
{
  const auto & vertices = face.vertices();
  const size_t ivertex = std::distance(vertices.begin(), std::find(
      vertices.begin(), vertices.end(), vertex));
  size_t iv1, iv2;
  if (ivertex < vertices.size() - 1)
      iv2 = ivertex + 1;
  else iv2 = 0;
  if (ivertex > 0)
      iv1 = ivertex - 1;
  else iv1 = vertices.size() - 1;
  
  return {{vertex, vertices[iv1]}, {vertex, vertices[iv2]}};
}

void IntegrationRuleVerticesAverage::build_tributary_shapes_cells_()
{
  const auto & cell = _element._parent_cell;
  const auto & vertices = cell.vertices();
  for (size_t ivertex = 0; ivertex < vertices.size(); ++ivertex)
  {
    std::vector<std::vector<size_t>> tributary_cell_faces;
    angem::PointSet<3, double> tributary_cell_vertices;
    const size_t vertex = vertices[ivertex];
    for (const auto face : cell.faces())
    {
      if (face->has_vertex(vertex))
      {
        const auto edges = get_edges_with_vertex(vertex, *face);
        build_tributary_cell_faces_(edges, tributary_cell_faces, tributary_cell_vertices);
      }
    }
    _tributary3d.emplace_back(tributary_cell_vertices.points, tributary_cell_faces);
  }
}

void IntegrationRuleVerticesAverage::
build_tributary_cell_faces_(const std::vector<mesh::Edge> & edges,
                            std::vector<std::vector<size_t>> & tributary_faces,
                            angem::PointSet<3,double> & tributary_vertices)
{
  const auto & grid = _element._parent_grid;
  const size_t v0 = tributary_vertices.insert(grid.vertex(edges[0].first));
  const Point p1 = grid.vertex(edges[0].first)  + grid.vertex(edges[0].second);
  const Point p2 = grid.vertex(edges[1].first)  + grid.vertex(edges[1].second);
  const size_t v1 = tributary_vertices.insert(p1);
  const size_t v2 = tributary_vertices.insert(p2);
  tributary_faces.push_back({v0, v1, v2});
  // cell center
  const Point c = _element._parent_cell.center();
  const size_t vc = tributary_vertices.insert(c);
  tributary_faces.push_back({vc, v1, v2});
}

void IntegrationRuleVerticesAverage::
build_tributary_shapes_face_(const size_t iface, const angem::Polygon<double> & face_poly)
{
  _tributary2d.emplace_back();
  auto & face_tributary_polygons = _tributary2d.back();
  const auto & points = face_poly.get_points();
  const size_t nv = points.size();
  const Point c = face_poly.center();
  for (size_t ivertex = 0; ivertex < nv; ++ivertex)
  {
    std::vector<Point> tributary_vertices(4);
    // first vertex
    tributary_vertices[0] = face_poly.get_points()[ivertex];
    // second vertex
    const size_t vertex1 = (ivertex < nv-1) ? ivertex + 1 : 0;
    tributary_vertices[1] = 0.5 * (tributary_vertices[0] + points[vertex1]);
    // third vertex
    const size_t vertex2 = (ivertex > 0) ? ivertex - 1 : nv - 1;
    tributary_vertices[2] = 0.5 * (tributary_vertices[0] + points[vertex2]);
    // face center
    tributary_vertices[3] = c;
    face_tributary_polygons.emplace_back(tributary_vertices);
  }
}

void IntegrationRuleVerticesAverage::compute_cell_fe_quantities_()
{
  const Point parent_center = _element._parent_cell.center();
  bool center_found = false;
  // integrate over regions
  const size_t n_parents = _element._basis_functions.size();
  std::vector<double> region_volumes( _tributary3d.size(), 0.0 );
  FeValues<angem::VTK_ID::TetrahedronID> fe_values;
  const auto & grid = _element._subgrid;
  auto & cell_data = _element._cell_data;
  for( auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell  )
  {
    const Point c = cell->center();
    fe_values.update(*cell);
    for (size_t region=0; region < _tributary3d.size(); ++region)  // tributary regions
    {
      if (_tributary3d[region].point_inside(c))
      {
        const std::vector<size_t> & cell_verts = cell->vertices();
        const size_t nv = cell_verts.size();
        region_volumes[region] += cell->volume();
        auto & data = cell_data.points[region];

        for (size_t parent_vertex=0; parent_vertex<n_parents; ++parent_vertex)
          for (size_t v=0; v<nv; ++v)
            for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
            {
              data.values[parent_vertex] += fe_values.value( v, q ) *
                  _element._basis_functions[parent_vertex][cell_verts[v]] *
                  fe_values.JxW(q);
              data.grads[parent_vertex] += fe_values.grad(v, q) *
                  _element._basis_functions[parent_vertex][cell_verts[v]] *
                  fe_values.JxW(q);
          }

        break;  // stop searching region
      }
    }
    const std::vector<size_t> & cell_verts = cell->vertices();
    const size_t nv = cell_verts.size();
    for (size_t parent_vertex=0; parent_vertex<n_parents; ++parent_vertex)
      for (size_t v=0; v<nv; ++v)
        for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
        {
          cell_data.center.values[parent_vertex] += fe_values.value(v, q) *
              _element._basis_functions[parent_vertex][cell_verts[v]] *
              fe_values.JxW(q);
          cell_data.center.grads[parent_vertex] += fe_values.grad(v, q) *
              _element._basis_functions[parent_vertex][cell_verts[v]] *
              fe_values.JxW(q);
        }
  }

  // normalize values and grads  by region volume
  for (size_t region=0; region < _tributary3d.size(); ++region)  // tributary regions
  {
    auto & data = cell_data.points[region];
    for (size_t parent_vertex=0; parent_vertex<n_parents; ++parent_vertex)
    {
      data.values[parent_vertex] /= region_volumes[region];
      data.grads[parent_vertex] /= region_volumes[region];
    }
    data.weight = region_volumes[region];
  }

  const double vol = _element._parent_cell.volume();
  auto & data = cell_data.center;
  for (size_t parent_vertex=0; parent_vertex<n_parents; ++parent_vertex)
  {
    data.values[parent_vertex] /= vol;
    data.grads[parent_vertex] /= vol;
  }
  data.weight = vol;
}

void IntegrationRuleVerticesAverage::compute_face_fe_quantities_(const size_t parent_face)
{
  const auto & grid = _element._subgrid;
  const std::vector<size_t> & face_indices = _element._face_domains[parent_face];
  FeValues<angem::VTK_ID::TriangleID> fe_values;
  const auto basis = grid
                     .face(face_indices.front())
                     .polygon()
                     .plane()
                     .get_basis();
  fe_values.set_basis(basis);
  const auto & regions = _tributary2d[parent_face];
  std::vector<double> region_areas( regions.size(), 0.0 );
  const size_t n_parent_vertices = _element._face_data[parent_face].center.values.size();
  // map parent face vertices to parent cell vertices
  // to get the index of the basis function
  std::vector<size_t> parent_vertices(n_parent_vertices);
  const auto face = _element._parent_cell.faces()[parent_face];
  const auto parent_cell_vertices = _element._parent_cell.vertices();
  for (size_t v=0; v<n_parent_vertices; ++v)
    parent_vertices[v] =
        std::distance(parent_cell_vertices.begin(),
                      std::find( parent_cell_vertices.begin(), parent_cell_vertices.end(),
                                 face->vertices()[v]));

  const auto & basis_functions = _element._basis_functions;

  size_t nhits = 0;
  for (const size_t iface : face_indices)
  {
    const mesh::Face & face = grid.face(iface);
    const Point c = face.center();
    fe_values.update(face);

    for (size_t region=0; region<regions.size(); ++region)  // tributary regions
      if (regions[region].point_inside(c))
      {
        nhits++;
        const std::vector<size_t> & face_verts = face.vertices();
        const size_t nv = face_verts.size();
        for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
          region_areas[region] += fe_values.JxW(q);

        auto & data = _element._face_data[parent_face].points[region];

        for (size_t parent_vertex=0; parent_vertex<n_parent_vertices; ++parent_vertex)
          for (size_t v=0; v<nv; ++v)
            for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
            {
              data.values[parent_vertex] += fe_values.value( v, q ) *
                  basis_functions[parent_vertices[parent_vertex]][face_verts[v]] *
                  fe_values.JxW(q);
              data.grads[parent_vertex] += fe_values.grad(v, q) *
                  basis_functions[parent_vertices[parent_vertex]][face_verts[v]];
                  fe_values.JxW(q);
          }

        break;  // stop searching region
      }

    const auto parent_cell_vertices = _element._parent_cell.vertices();
    const std::vector<size_t> & face_verts = face.vertices();
    const size_t nv = face_verts.size();
    auto & data = _element._face_data[parent_face].center;
    for (size_t parent_vertex=0; parent_vertex<n_parent_vertices; ++parent_vertex)
      for (size_t v=0; v<nv; ++v)
        for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
        {
          data.values[parent_vertex] += fe_values.value(v, q) *
                      basis_functions[parent_vertices[parent_vertex]][face_verts[v]] *
                      fe_values.JxW(q);
          data.grads[parent_vertex] += fe_values.grad(v, q) *
                      basis_functions[parent_vertices[parent_vertex]][face_verts[v]] *
                      fe_values.JxW(q);
        }
    data.weight += face.area();
  }

  if (nhits != face_indices.size())
    throw std::runtime_error("small nhits");

  // normalize values and grads  by region volume
  double sum_weights = 0;
  for (size_t region=0; region<regions.size(); ++region)  // tributary regions
  {
    auto & data = _element._face_data[parent_face].points[region];
    for (size_t parent_vertex=0; parent_vertex<n_parent_vertices; ++parent_vertex)
    {
      if (std::isnan(data.weight) || std::fabs(region_areas[region]) < 1e-10)
      {
        std::cout << "\n" << std::endl << std::flush;
        auto f = _element._parent_cell.faces()[parent_face];
        {
          for (auto v : f->vertices())
            std::cout << v << " ";
          std::cout << std::endl;
        }
        std::cout << "region = " << region << std::endl;
        std::cout << "regions.size() = " << regions.size() << std::endl;
        for (size_t r=0; r<regions.size(); ++r)
          std::cout << "area = " << region_areas[r] << std::endl << std::flush;

        _element.save_shape_functions("output/mess-" + std::to_string(_element._parent_cell.index()) + ".vtk");
        std::cout << "_element_parent_cell.index() = " << _element._parent_cell.index() << std::endl;
        std::cout << "region area = " << region_areas[region] << std::endl << std::flush;
        throw std::runtime_error("bad face integration weight");
      }
      data.values[parent_vertex] /= region_areas[region];
      data.grads[parent_vertex] /= region_areas[region];
    }
    data.weight = region_areas[region];
    sum_weights += region_areas[region];
  }
  if (std::fabs(sum_weights - _element._parent_cell.faces()[parent_face]->area()) > 1e-10)
    throw std::runtime_error("weights and face area don't add up");

  auto & data = _element._face_data[parent_face].center;
  for (size_t parent_vertex=0; parent_vertex<n_parent_vertices; ++parent_vertex)
  {
    data.values[parent_vertex] /= data.weight;
    data.grads[parent_vertex] /= data.weight;
  }
}



}  // end namespace discretization

#endif
