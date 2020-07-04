#ifdef WITH_EIGEN
#include "IntegrationRuleFractureFull.hpp"
#include "../FeValues.hpp"

namespace discretization {

IntegrationRuleFractureFull::IntegrationRuleFractureFull(PolyhedralElementBase & element)
    : _element(element)
{
  if (_element._face_domains.empty())
    _element._face_domains = _element.create_face_domains_(/*sort_faces=*/true);

  setup_storage_();
  const size_t n_faces = _element._face_domains.size();
  for (size_t iface=0; iface < n_faces; ++iface)
    compute_face_fe_quantities_(iface);
}

void IntegrationRuleFractureFull::setup_storage_()
{
  auto & data = _element._face_fracture_data;
  data.resize(_element._parent_cell.faces().size());
  const size_t n_parent_vertices = _element._parent_cell.vertices().size();
  for (size_t iface=0; iface<data.size(); ++iface)
  {
    const size_t nq = _element._face_domains[iface].size();
    data[iface].points.resize(nq);
    for (size_t q = 0; q < nq; ++q)
    {
      data[iface].points[q].values.resize( n_parent_vertices );
      data[iface].points[q].grads.resize( n_parent_vertices );
    }
  }
}

namespace {  // anonymous NS since the function is defined in another cpp file

void get_face_integration_points(FeValues<angem::TriangleID> & fe_values,
                                 std::vector<Point> & integration_points)
{
  const std::vector<Point> master_qpoints = fe_values.get_master_integration_points();
  if (integration_points.size() != master_qpoints.size())
    integration_points.resize( master_qpoints.size() );
  Point zero_point = {0,0,0};
  std::fill(integration_points.begin(), integration_points.end(), zero_point);

  const size_t nv = ElementTraits<angem::TriangleID>::n_vertices;  // number of vertices in triangle
  for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
    for (size_t v = 0; v < nv; ++v)
      integration_points[q] += fe_values.value(v, q) * master_qpoints[q];
}

}

void IntegrationRuleFractureFull::compute_face_fe_quantities_(const size_t parent_face)
{
  const auto & grid = _element._element_grid;
  const std::vector<size_t> & face_indices = _element._face_domains[parent_face];
  FeValues<angem::VTK_ID::TetrahedronID> fe_cell_values;
  FeValues<angem::VTK_ID::TriangleID> fe_face_values;
  const auto basis = grid
                     .face(face_indices.front())
                     .polygon()
                     .plane()
                     .get_basis();

  fe_face_values.set_basis(basis);

  std::vector<Point> local_integration_points;
  const size_t n_parent_vertices = _element._parent_cell.vertices().size();
  const size_t nv = ElementTraits<angem::TetrahedronID>::n_vertices;
  const auto & basis_functions = _element._basis_functions;

  for (size_t iface = 0; iface < face_indices.size(); ++iface)
  {
    const mesh::Face & face = grid.face(face_indices[iface]);
    fe_face_values.update(face);
    get_face_integration_points(fe_face_values, local_integration_points);
    const mesh::Cell & cell = *face.neighbors()[0];  // face only has one neighbor
    fe_cell_values.update(cell, local_integration_points);

    const std::vector<size_t> & cell_verts = cell.vertices();
    auto & data = _element._face_fracture_data[parent_face].points[iface];

    for (size_t parent_vertex = 0; parent_vertex < n_parent_vertices; ++parent_vertex)
      for (size_t v = 0; v < nv; ++v)
      {
        data.values[parent_vertex] += fe_cell_values.value(v, 0) *
                                     basis_functions[parent_vertex][cell_verts[v]];
        data.grads[parent_vertex] += fe_cell_values.grad(v, 0) *
                                     basis_functions[parent_vertex][cell_verts[v]];
      }
    data.weight = face.area();
  }
}


}  // end namespace discretization


#endif