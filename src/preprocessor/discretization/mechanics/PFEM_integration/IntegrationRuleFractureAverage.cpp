#ifdef WITH_EIGEN
#include "IntegrationRuleFractureAverage.hpp"
#include "../FeValues.hpp"

namespace discretization {

using Point = angem::Point<3,double>;

IntegrationRuleFractureAverage::
IntegrationRuleFractureAverage(PolyhedralElementBase & element,
                               const std::vector<std::vector<angem::Polygon<double>>> & tributary_2d,
                               const size_t parent_face)
    : _element(element), _tributary_2d(tributary_2d), _parent_face(parent_face)
{
  if (_element._face_domains.empty())
    _element._face_domains = _element.create_face_domains_();
}

void IntegrationRuleFractureAverage::setup_storage_(FiniteElementData & data) const
{
  // auto & data = _element._face_fracture_data;
  // data.resize(_element._parent_cell.faces().size());
  const size_t n_parent_vertices = _element._parent_cell.vertices().size();
  // for (size_t iface=0; iface<data.size(); ++iface)
  const size_t nq = _tributary_2d[_parent_face].size();
  data.points.resize(nq);
  for (size_t q = 0; q < nq; ++q)
  {
    data.points[q].values.resize( n_parent_vertices );
    data.points[q].grads.resize( n_parent_vertices );
  }
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
  const auto & grid = _element._element_grid;
  const std::vector<size_t> & face_indices = _element._face_domains[_parent_face];
  FeValues<angem::VTK_ID::TetrahedronID> fe_cell_values;
  FeValues<angem::VTK_ID::TriangleID> fe_face_values;
  const auto basis = grid.face(face_indices.front()).polygon().plane().get_basis();
  fe_face_values.set_basis(basis);

  const auto & regions = _tributary_2d[_parent_face];
  std::vector<double> region_areas( regions.size(), 0.0 );
  const size_t n_parent_vertices = _element._parent_cell.vertices().size();
  const auto & basis_functions = _element._basis_functions;
  std::vector<Point> local_integration_points;
  const size_t nv = ElementTraits<angem::TetrahedronID>::n_vertices;


  for (const size_t iface : face_indices)
  {
    const mesh::Face & face = grid.face(iface);
    const Point c = face.center();

    for (size_t region=0; region<regions.size(); ++region)  // tributary regions
      if (regions[region].point_inside(c))
      {
        const mesh::Cell & cell = *face.neighbors()[0];  // face only has one neighbor
        const std::vector<size_t> & cell_verts = cell.vertices();
       
        fe_face_values.update(face);
        get_face_integration_points(fe_face_values, local_integration_points);
        fe_cell_values.update(cell, local_integration_points);
        // fe_cell_values.update(cell);

        for (size_t q = 0; q < fe_face_values.n_integration_points(); ++q)
          region_areas[region] += fe_face_values.JxW(q);

        auto & data = face_data.points[region];

        for (size_t parent_vertex = 0; parent_vertex < n_parent_vertices; ++parent_vertex)
          for (size_t v = 0; v < nv; ++v)
            for (size_t q = 0; q < fe_cell_values.n_integration_points(); ++q)
            {
              data.values[parent_vertex] += fe_cell_values.value(v, q) *
                                            basis_functions[parent_vertex][cell_verts[v]] *
                                            fe_face_values.JxW(q);
              data.grads[parent_vertex] += fe_cell_values.grad(v, q) *
                                           basis_functions[parent_vertex][cell_verts[v]] *
                                           fe_face_values.JxW(q);
            }

        break;  // stop searching region
      }
  }

  for (size_t region=0; region<regions.size(); ++region)  // tributary regions
  {
    auto & data = face_data.points[region];
    for (size_t parent_vertex=0; parent_vertex<n_parent_vertices; ++parent_vertex)
    {
      data.values[parent_vertex] /= region_areas[region];
      data.grads[parent_vertex] /= region_areas[region];
    }
  }
  return face_data;
}


}  // end namespace discretization

#endif
