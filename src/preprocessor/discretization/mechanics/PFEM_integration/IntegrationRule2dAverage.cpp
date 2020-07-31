#include "IntegrationRule2dAverage.hpp"
#include "TributaryRegion2dBase.hpp"
#include "../FeValues.hpp"

#ifdef WITH_EIGEN

namespace discretization {

IntegrationRule2dAverage::IntegrationRule2dAverage(PolyhedralElementBase & element,
                           const TributaryRegion2dBase & tributary,
                           const size_t parent_face,
                           const angem::Basis<3, double> & basis)
    : IntegrationRule2dBase(element, parent_face, basis), _tributary(tributary)
{
  assert( tributary.size() > 0 );
}

void IntegrationRule2dAverage::setup_storage_(FiniteElementData & data) const
{
  const size_t n_vertices = compute_n_parent_vertices_();
  data.points.resize( _tributary.size() );
  for (size_t q=0; q<data.points.size(); ++q)
  {
    data.points[q].values.resize(n_vertices);
    data.points[q].grads.resize(n_vertices);
  }
  data.center.values.resize(n_vertices);
  data.center.grads.resize(n_vertices);
}

void IntegrationRule2dAverage::compute_fe_values_(const std::vector<size_t> &faces, FEPointData &dst) const
{
  const auto & grid = _element._subgrid;
  FeValues<angem::VTK_ID::TriangleID> fe_values;
  fe_values.set_basis(_basis);
  double sum_areas = 0;
  const auto & basis_functions = _element._basis_functions;
  for (const size_t iface : faces)
  {
    const mesh::Face & face = grid.face(iface);
    fe_values.update(face);
    const std::vector<size_t> & face_verts = face.vertices();

    for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
      sum_areas += fe_values.JxW(q);

    for (size_t parent_vertex=0; parent_vertex < _parent_vertices.size(); ++parent_vertex)
      for (size_t v=0; v<face_verts.size(); ++v)
        for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
        {
          dst.values[parent_vertex] += fe_values.value( v, q ) *
              basis_functions[_parent_vertices[parent_vertex]][face_verts[v]] *
              fe_values.JxW(q);
          dst.grads[parent_vertex] += fe_values.grad(v, q) *
              basis_functions[_parent_vertices[parent_vertex]][face_verts[v]] *
              fe_values.JxW(q);
        }
  }

  std::transform(dst.values.begin(), dst.values.end(), dst.values.begin(),
                 [sum_areas](const double value) {return value / sum_areas;} );
  std::transform(dst.grads.begin(), dst.grads.end(), dst.grads.begin(),
                 [sum_areas](const Point & value) {return value / sum_areas;} );
}

FiniteElementData IntegrationRule2dAverage::get() const
{
  FiniteElementData face_data;
  setup_storage_(face_data);

  for (size_t region = 0; region < _tributary.size(); ++region)
  {
    compute_fe_values_(_tributary.get(region), face_data.points[region]);
    face_data.points[region].weight = _tributary.area(region);
  }

  compute_fe_values_(_tributary.get_center(), face_data.center);
  face_data.center.weight = _tributary.area_total();
  return face_data;
}

}  // end namespace discretization

#endif
