#include "IntegrationRuleFacesAverage.hpp"
#include "../FeValues.hpp"

namespace discretization {

using Point = angem::Point<3,double>;

IntegrationRuleFacesAverage::IntegrationRuleFacesAverage(PolyhedralElementBase & element)
    : _element(element)
{
  if (_element._face_domains.empty())
    _element._face_domains = _element.create_face_domains_();

  build_tributary_shapes_cells_();

  const auto face_polygons = _element._parent_cell.polyhedron()->get_face_polygons();
  const size_t n_faces = face_polygons.size();
  _face_triangles.resize(n_faces);
  for (size_t iface=0; iface < n_faces; ++iface)
    build_tributary_shapes_face_(iface, face_polygons[iface]);

  setup_storage_();
  compute_cell_fe_quantities_();
  for (size_t iface=0; iface < n_faces; ++iface)
    compute_face_fe_quantities_(iface);
}

void IntegrationRuleFacesAverage::build_tributary_shapes_cells_()
{
  /* Split a parent cell into tributary regions (pyramids) */
  const auto polyhedron = _element._parent_cell.polyhedron();
  std::vector<Point> vertices = polyhedron->get_points();
  vertices.push_back(_element._parent_cell.center());
  _element._cell_gauss_points.clear();
  for (const auto & face : polyhedron->get_faces())
  {
    _pyramids.push_back(create_pyramid_(face, vertices));
    _element._cell_gauss_points.push_back( _pyramids.back().center() );
  }
}

angem::Polyhedron<double>
IntegrationRuleFacesAverage::create_pyramid_(const std::vector<size_t> & face,
                                             const std::vector<Point> & vertices) const
{
  const size_t vertex_center = vertices.size() - 1;  // HACK: I just pushed it to this array
  const auto c = _element._parent_cell.center();
  std::vector<std::vector<size_t>> pyramid_faces;
  for (size_t iv=0; iv<face.size(); ++iv)
  {
    size_t v1, v2;
    v1 = face[iv];
    if ( iv < face.size() - 1 )
      v2 = face[iv+1];
    else
      v2 = face[0];
    pyramid_faces.push_back( {v1, v2, vertex_center} );
  }
  pyramid_faces.push_back( face );  // base

  angem::Polyhedron<double> pyramid(vertices, pyramid_faces);
  return pyramid;
}

void IntegrationRuleFacesAverage::compute_cell_fe_quantities_()
{
  const Point parent_center = _element._parent_cell.center();
  bool center_found = false;
  // integrate over regions
  const size_t n_parents = _element._basis_functions.size();
  std::vector<double> region_volumes( _pyramids.size(), 0.0 );
  FeValues<angem::VTK_ID::TetrahedronID> fe_values;
  const auto & grid = _element._element_grid;
  auto & cell_data = _element._cell_data;
  for( auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell  )
  {
    const Point c = cell->center();
    fe_values.update(*cell);
    for (size_t region=0; region<_pyramids.size(); ++region)  // tributary regions
    {
      if (_pyramids[region].point_inside(c))
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
  for (size_t region=0; region<_pyramids.size(); ++region)  // tributary regions
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

void IntegrationRuleFacesAverage::setup_storage_()
{
  auto & cell_data = _element._cell_data;
  const auto & basis_functions = _element._basis_functions;
  cell_data.points.resize( _pyramids.size() );
  for (auto & point : cell_data.points)
  {
    point.values.resize( basis_functions.size(), 0 );
    point.grads.resize( basis_functions.size() );
  }

  cell_data.center.values.resize( basis_functions.size(), 0 );
  cell_data.center.grads.resize( basis_functions.size() );

  const std::vector<const mesh::Face*> parent_faces = _element._parent_cell.faces();
  _element._face_data.resize( _face_triangles.size() );
  for (size_t parent_face = 0; parent_face < _face_triangles.size(); ++parent_face)
  {
    auto & data = _element._face_data[parent_face];
    const size_t n_vertices = parent_faces[parent_face]->vertices().size();
    data.points.resize( _face_triangles[parent_face].size() );
    for (size_t q=0; q<data.points.size(); ++q)
    {
      data.points[q].values.resize(n_vertices);
      data.points[q].grads.resize(n_vertices);
    }
    data.center.values.resize(n_vertices);
    data.center.grads.resize(n_vertices);
  }
}

void IntegrationRuleFacesAverage::build_tributary_shapes_face_(const size_t iface, const angem::Polygon<double> & face_poly)
{
  const angem::Point<3,double> center = face_poly.center();
  const std::vector<angem::Point<3,double>> & vertices = face_poly.get_points();
  std::vector<angem::Point<3,double>> result;
  for (size_t i=0; i<vertices.size(); ++i)
  {
    size_t v1 = i, v2;
    if ( i <  vertices.size() - 1) v2 = i + 1;
    else                           v2 = 0;
    std::vector<angem::Point<3,double>> triangle_vertices = {vertices[v1], vertices[v2], center};
    angem::Polygon<double> triangle(triangle_vertices);
    _face_triangles[iface].push_back(triangle);
  }
}

void IntegrationRuleFacesAverage::compute_face_fe_quantities_(const size_t parent_face)
{
  const auto & grid = _element._element_grid;
  const std::vector<size_t> & face_indices = _element._face_domains[parent_face];
  FeValues<angem::VTK_ID::TriangleID> fe_values;
  const auto basis = grid
                     .face(face_indices.front())
                     .polygon()
                     .plane()
                     .get_basis();
  fe_values.set_basis(basis);
  const auto & regions = _face_triangles[parent_face];
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
