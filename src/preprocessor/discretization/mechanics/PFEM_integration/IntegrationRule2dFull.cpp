#include "IntegrationRule2dFull.hpp"
#include "../FeValues.hpp"

#ifdef WITH_EIGEN

namespace discretization {

IntegrationRule2dFull::IntegrationRule2dFull(PolyhedralElementBase & element, const size_t parent_face,
                        const angem::Basis<3, double> & basis)
    : _element(element), _parent_face(parent_face), _basis(basis)
{
  if (_element._face_domains.empty())
    _element._face_domains = _element.create_face_domains_();
}

void IntegrationRule2dFull::setup_storage_(FiniteElementData & data) const
{
  const std::vector<const mesh::Face*> parent_faces = _element._parent_cell.faces();
    data.points.resize(_element._face_domains[_parent_face].size());
    const size_t n_vertices = parent_faces[_parent_face]->vertices().size();
    for (size_t q = 0; q < data.points.size(); ++q)
    {
      data.points[q].values.resize(n_vertices);
      data.points[q].grads.resize(n_vertices);
    }
    data.center.values.resize(n_vertices);
    data.center.grads.resize(n_vertices);
}

FiniteElementData IntegrationRule2dFull::get() const
{
  FiniteElementData face_data;
  setup_storage_(face_data);
  const auto & grid = _element._subgrid;
  const size_t ipf = _parent_face;
  const std::vector<size_t> & face_indices = _element._face_domains[_parent_face];
  FeValues<angem::VTK_ID::TriangleID> fe_values;
  fe_values.set_basis(_basis);
  const size_t n_parent_vertices = face_data.center.values.size();
  const auto * parent_face = _element._parent_cell.faces()[_parent_face];
  const Point parent_center = parent_face->center();
  const auto parent_cell_vertices = _element._parent_cell.vertices();
  std::vector<size_t> parent_vertices(n_parent_vertices);
  double min_dist = std::numeric_limits<double>::max();
  for (size_t v=0; v<n_parent_vertices; ++v)
    parent_vertices[v] = std::distance(parent_cell_vertices.begin(),
                                std::find( parent_cell_vertices.begin(), parent_cell_vertices.end(),
                                 parent_face->vertices()[v]));
  const auto & basis_functions = _element._basis_functions;

  auto & data_center = face_data.center;
  for (size_t iface = 0; iface < face_indices.size(); ++iface)
  {
    const size_t face_index = face_indices[iface];
    const auto & face = grid.face(face_index);
    const auto & face_verts = face.vertices();
    const size_t nv = face_verts.size();
    auto & data = face_data.points[iface];
    fe_values.update(face);
    for (size_t parent_vertex=0; parent_vertex < n_parent_vertices; ++parent_vertex)
      for (size_t v=0; v<nv; ++v)
      {
        data.values[parent_vertex] += fe_values.value( v, 0 ) *
                                      basis_functions[parent_vertices[parent_vertex]][face_verts[v]];
        data.grads[parent_vertex] += fe_values.grad(v, 0) *
                                     basis_functions[parent_vertices[parent_vertex]][face_verts[v]];
        data_center.values[parent_vertex] += fe_values.value( v, 0 ) *
                                             basis_functions[parent_vertices[parent_vertex]][face_verts[v]] *
                                             fe_values.JxW(0);
        data_center.grads[parent_vertex] += fe_values.grad( v, 0 ) *
                                            basis_functions[parent_vertices[parent_vertex]][face_verts[v]] *
                                            fe_values.JxW(0);
      }
    data.weight = face.area();
  }
  const double parent_face_area = parent_face->area();
  face_data.center.weight = parent_face_area;
  for (size_t parent_vertex=0; parent_vertex < n_parent_vertices; ++parent_vertex)
  {
    face_data.center.values[parent_vertex] /= parent_face_area;
    face_data.center.grads[parent_vertex] /= parent_face_area;
  }

  return face_data;
}

}  // end namespace discretization


#endif
