#include "IntegrationRule2dAverage.hpp"
#include "../FeValues.hpp"

#ifdef WITH_EIGEN

namespace discretization {

IntegrationRule2dAverage::IntegrationRule2dAverage(PolyhedralElementBase & element,
                           const std::vector<std::vector<angem::Polygon<double>>> & tributary_2d,
                           const size_t parent_face,
                           const angem::Basis<3, double> & basis)
    : _element(element), _tributary_2d(tributary_2d), _parent_face(parent_face), _basis(basis)
{
  if (_element._face_domains.empty())
    _element._face_domains = _element.create_face_domains_();
}

void IntegrationRule2dAverage::setup_storage_(FiniteElementData & data) const
{
  const std::vector<const mesh::Face*> parent_faces = _element._parent_cell.faces();
  const size_t n_vertices = parent_faces[_parent_face]->vertices().size();
  data.points.resize( _tributary_2d[_parent_face].size() );
  for (size_t q=0; q<data.points.size(); ++q)
  {
    data.points[q].values.resize(n_vertices);
    data.points[q].grads.resize(n_vertices);
  }
  data.center.values.resize(n_vertices);
  data.center.grads.resize(n_vertices);
}

FiniteElementData IntegrationRule2dAverage::get() const
{
  FiniteElementData face_data;
  setup_storage_(face_data);
  const auto & grid = _element._subgrid;
  const std::vector<size_t> & face_indices = _element._face_domains[_parent_face];
  FeValues<angem::VTK_ID::TriangleID> fe_values;
  // const auto basis = grid.face(face_indices.front()).polygon().plane().get_basis();
  fe_values.set_basis(_basis);
  const auto & regions = _tributary_2d[_parent_face];
  std::vector<double> region_areas( regions.size(), 0.0 );
  const size_t n_parent_vertices = face_data.center.values.size();

  // map parent face vertices to parent cell vertices
  // to get the index of the basis function
  std::vector<size_t> parent_vertices(n_parent_vertices);
  const auto face = _element._parent_cell.faces()[_parent_face];
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

        auto & data = face_data.points[region];

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
    auto & data = face_data.center;
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
    auto & data = face_data.points[region];
    for (size_t parent_vertex=0; parent_vertex<n_parent_vertices; ++parent_vertex)
    {
      if (std::isnan(data.weight) || std::fabs(region_areas[region]) < 1e-10)
      {
        std::cout << "\n" << std::endl << std::flush;
        auto f = _element._parent_cell.faces()[_parent_face];
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
  if (std::fabs(sum_weights - _element._parent_cell.faces()[_parent_face]->area()) > 1e-10)
  {
    std::cout << std::endl;
    std::cout << "geom area = " << _element._parent_cell.faces()[_parent_face]->area() << std::endl;
    for (auto w : region_areas)
      std::cout << w << " ";
    std::cout << " (sum) " << sum_weights << std::endl;
    throw std::runtime_error("weights and face area don't add up");
  }

  auto & data = face_data.center;
  for (size_t parent_vertex=0; parent_vertex<n_parent_vertices; ++parent_vertex)
  {
    data.values[parent_vertex] /= data.weight;
    data.grads[parent_vertex] /= data.weight;
  }
  return face_data;
}

}  // end namespace discretization

#endif
