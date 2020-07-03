#include "IntegrationRule2dFull.hpp"
#include "../FeValues.hpp"

#ifdef WITH_EIGEN

namespace discretization {

IntegrationRule2dFull::IntegrationRule2dFull(PolyhedralElementBase & element)
    : _element(element)
{
  if (_element._face_domains.empty())
    _element._face_domains = _element.create_face_domains_();
  setup_storage_();
  const size_t n_faces = _element._face_domains.size();
  for (size_t iface=0; iface < n_faces; ++iface)
    compute_face_fe_quantities_(iface);
}

void IntegrationRule2dFull::setup_storage_()
{
  const std::vector<const mesh::Face*> parent_faces = _element._parent_cell.faces();
  // const auto & basis_functions = _element._basis_functions;
  _element._face_data.resize(_element._face_domains.size());
  for (size_t parent_face = 0; parent_face < _element._face_domains.size(); ++parent_face)
  {
    auto & data = _element._face_data[parent_face];
    data.points.resize(_element._face_domains[parent_face].size());
    const size_t n_vertices = parent_faces[parent_face]->vertices().size();
    for (size_t q = 0; q < data.points.size(); ++q)
    {
      data.points[q].values.resize(n_vertices);
      data.points[q].grads.resize(n_vertices);
    }
    data.center.values.resize(n_vertices);
    data.center.grads.resize(n_vertices);
  }
}

void IntegrationRule2dFull::compute_face_fe_quantities_(const size_t ipf)
{
  const auto & grid = _element._element_grid;
  const std::vector<size_t> & face_indices = _element._face_domains[ipf];
  FeValues<angem::VTK_ID::TriangleID> fe_values;
  const auto basis = grid.face(face_indices.front()).polygon().plane().get_basis();
  fe_values.set_basis(basis);
  const size_t n_parent_vertices = _element._face_data[ipf].center.values.size();
  const auto * parent_face = _element._parent_cell.faces()[ipf];
  const Point parent_center = parent_face->center();
  const auto parent_cell_vertices = _element._parent_cell.vertices();
  std::vector<size_t> parent_vertices(n_parent_vertices);
  double min_dist = std::numeric_limits<double>::max();
  for (size_t v=0; v<n_parent_vertices; ++v)
    parent_vertices[v] = std::distance(parent_cell_vertices.begin(),
                                std::find( parent_cell_vertices.begin(), parent_cell_vertices.end(),
                                 parent_face->vertices()[v]));
  const auto & basis_functions = _element._basis_functions;

  for (size_t iface = 0; iface < face_indices.size(); ++iface)
  {
    const size_t face_index = face_indices[iface];
    const auto & face = grid.face(face_index);
    const auto & face_verts = face.vertices();
    const size_t nv = face_verts.size();
    auto & data = _element._face_data[ipf].points[iface];
    fe_values.update(face);
    for (size_t parent_vertex=0; parent_vertex < n_parent_vertices; ++parent_vertex)
      for (size_t v=0; v<nv; ++v)
      {
        data.values[parent_vertex] += fe_values.value( v, 0 ) *
            basis_functions[parent_vertex][face_verts[v]];
        data.grads[parent_vertex] += fe_values.grad(v, 0) *
            basis_functions[parent_vertex][face_verts[v]];
        _element._face_data[ipf].center.values[parent_vertex] +=
            fe_values.value( v, 0 ) *
            basis_functions[parent_vertex][face_verts[v]] *
            fe_values.JxW(0);
        _element._face_data[ipf].center.grads[parent_vertex] +=
            fe_values.grad( v, 0 ) *
            basis_functions[parent_vertex][face_verts[v]] *
            fe_values.JxW(0);
      }
    data.weight = face.area();
  }
  const double parent_face_area = parent_face->area();
  _element._face_data[ipf].center.weight = parent_face_area;
  for (size_t parent_vertex=0; parent_vertex < n_parent_vertices; ++parent_vertex)
  {
    _element._face_data[ipf].center.values[parent_vertex] /= parent_face_area;
    _element._face_data[ipf].center.grads[parent_vertex] /= parent_face_area;
  }
}

}  // end namespace discretization


#endif
