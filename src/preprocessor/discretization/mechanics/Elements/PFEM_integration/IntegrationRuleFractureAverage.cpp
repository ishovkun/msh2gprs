#ifdef WITH_EIGEN
#include "IntegrationRuleFractureAverage.hpp"
#include "TributaryRegion2dBase.hpp"
#include "../FeValues.hpp"

namespace discretization {

using Point = angem::Point<3,double>;

IntegrationRuleFractureAverage::
IntegrationRuleFractureAverage(PolyhedralElementBase & element,
                               const TributaryRegion2dBase & tributary_2d,
                               const size_t parent_face,
                               const angem::Basis<3, double> & basis)
    : _element(element), _tributary(tributary_2d), _parent_face(parent_face), _basis(basis)
{
  if (_element._face_domains.empty())
    _element._face_domains = _element.create_face_domains_();
}

void IntegrationRuleFractureAverage::setup_storage_(FiniteElementData & data) const
{
  const size_t n_parent_vertices = _element._parent_cell.vertices().size();
  data.points.resize(_tributary.size());
  for (size_t q = 0; q < _tributary.size(); ++q)
  {
    data.points[q].values.resize( n_parent_vertices, 0.0 );
    data.points[q].grads.resize( n_parent_vertices, {0,0,0} );
  }
  data.center.values.resize(n_parent_vertices, 0.0);
  data.center.grads.resize(n_parent_vertices, {0,0,0});
}

void get_face_integration_points(FeValues<angem::TriangleID> & fe_values,
                                 std::vector<Point> & integration_points)
{
  const std::vector<Point> master_qpoints = fe_values.get_master_integration_points();
  if (integration_points.size() != master_qpoints.size())
    integration_points.resize( master_qpoints.size() );
  Point zero_point = {0,0,0};
  std::fill(integration_points.begin(), integration_points.end(), zero_point);

  // const size_t nv = N_ELEMENT_VERTICES<angem::TriangleID>;  // number of vertices in triangle
  const size_t nv = ElementTraits<angem::TriangleID>::n_vertices;  // number of vertices in triangle
  for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
    for (size_t v = 0; v < nv; ++v)
      integration_points[q] += fe_values.value(v, q) * master_qpoints[q];
}

FiniteElementData IntegrationRuleFractureAverage::get() const
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

void IntegrationRuleFractureAverage::compute_fe_values_(const std::vector<size_t> &faces, FEPointData &dst) const
{
  const auto & grid = _element._subgrid;
  FeValues<angem::VTK_ID::TetrahedronID> fe_cell_values;
  FeValues<angem::VTK_ID::TriangleID> fe_face_values;
  fe_face_values.set_basis(_basis);

  const size_t n_parent_vertices = _element._parent_cell.vertices().size();
  const auto & basis_functions = _element._basis_functions;
  std::vector<Point> local_integration_points;
  const size_t nv = ElementTraits<angem::TetrahedronID>::n_vertices;
  double sum_areas = 0;

  for (const size_t iface : faces)
  {
    const mesh::Face & face = grid.face(iface);
    const mesh::Cell & cell = *face.neighbors()[0];  // face only has one neighbor
    const std::vector<size_t> & cell_verts = cell.vertices();
       
    fe_face_values.update(face);
    get_face_integration_points(fe_face_values, local_integration_points);
    fe_cell_values.update(cell, local_integration_points);

    for (size_t q = 0; q < fe_face_values.n_integration_points(); ++q)
      sum_areas += fe_face_values.JxW(q);

    for (size_t parent_vertex = 0; parent_vertex < n_parent_vertices; ++parent_vertex)
      for (size_t v = 0; v < nv; ++v)
        for (size_t q = 0; q < fe_cell_values.n_integration_points(); ++q)
        {
          dst.values[parent_vertex] += fe_cell_values.value(v, q) *
              basis_functions[parent_vertex][cell_verts[v]] *
              fe_face_values.JxW(q);
          dst.grads[parent_vertex] += fe_cell_values.grad(v, q) *
              basis_functions[parent_vertex][cell_verts[v]] *
              fe_face_values.JxW(q);
        }

  }

  std::transform(dst.values.begin(), dst.values.end(), dst.values.begin(),
                 [sum_areas](const double val) {return val / sum_areas;});
  std::transform(dst.grads.begin(), dst.grads.end(), dst.grads.begin(),
                 [sum_areas](const Point & val) {return val / sum_areas;});
}


}  // end namespace discretization

#endif
