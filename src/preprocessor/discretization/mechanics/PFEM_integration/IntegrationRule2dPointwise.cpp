#include "IntegrationRule2dPointwise.hpp"
#include "../FeValues.hpp"

#ifdef WITH_EIGEN

namespace discretization {

IntegrationRule2dPointwise::IntegrationRule2dPointwise(PolyhedralElementBase & element,
                           const std::vector<std::vector<angem::Polygon<double>>> & tributary_2d,
                           const size_t parent_face)
    : _element(element), _tributary_2d(tributary_2d), _parent_face(parent_face)
{
  if (_element._face_domains.empty())
    _element._face_domains = _element.create_face_domains_();
}

void IntegrationRule2dPointwise::setup_storage_(FiniteElementData & data) const
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

FiniteElementData IntegrationRule2dPointwise::get() const
{
  FiniteElementData face_data;
  setup_storage_(face_data);
  const auto & grid = _element._subgrid;
  const std::vector<size_t> & face_indices = _element._face_domains[_parent_face];
  FeValues<angem::VTK_ID::TriangleID> fe_values;
  const auto basis = grid.face(face_indices.front()).polygon().plane().get_basis();
  fe_values.set_basis(basis);
  const auto & regions = _tributary_2d[_parent_face];
  std::vector<int> region_found(regions.size(), 0);
  const size_t n_parent_vertices = face_data.center.values.size();
  // map parent face vertices to parent cell vertices
  // to get the index of the basis function
  std::vector<size_t> parent_vertices(n_parent_vertices);
  const auto * face = _element._parent_cell.faces()[_parent_face];
  const Point parent_center = face->center();

  const auto parent_cell_vertices = _element._parent_cell.vertices();
  for (size_t v=0; v<n_parent_vertices; ++v)
    parent_vertices[v] =
        std::distance(parent_cell_vertices.begin(),
                      std::find( parent_cell_vertices.begin(), parent_cell_vertices.end(),
                                 face->vertices()[v]));
  const auto & basis_functions = _element._basis_functions;

  bool center_found = false;
  for (const size_t iface : face_indices)
  {
    const mesh::Face & face = grid.face(iface);
    const auto face_poly = face.polygon();

    const std::vector<size_t> & face_verts = face.vertices();
    const size_t nv = face_verts.size();

    for (size_t region=0; region<regions.size(); ++region)  // tributary regions
      if (!region_found[region])
      if ( face_poly.point_inside(regions[region].center()) )
      {
        fe_values.update(face);
        region_found[region] = 1;
        auto & data = face_data.points[region];
        for (size_t parent_vertex=0; parent_vertex<n_parent_vertices; ++parent_vertex)
          for (size_t v=0; v<nv; ++v)
            {
              data.values[parent_vertex] += fe_values.value( v, 0 ) *
                  basis_functions[parent_vertex][face_verts[v]];
              data.grads[parent_vertex] += fe_values.grad(v, 0) *
                  basis_functions[parent_vertex][face_verts[v]];
            }
      }

    if (!center_found)
      if (face_poly.point_inside(parent_center))
      {
        center_found = true;
        fe_values.update(face, {parent_center});
        auto & data = face_data.center;
        const size_t q = 0;
        for (size_t parent_vertex=0; parent_vertex<n_parent_vertices; ++parent_vertex)
          for (size_t v=0; v<nv; ++v)
            {
              data.values[parent_vertex] += fe_values.value( v, q ) *
                  basis_functions[parent_vertex][face_verts[v]];
              data.grads[parent_vertex] += fe_values.grad(v, q) *
                  basis_functions[parent_vertex][face_verts[v]];
            }
      }
  }

  // normalize values and grads  by region volume
  for (size_t region=0; region<regions.size(); ++region)  // tributary regions
    face_data.points[region].weight = regions[region].area();
  face_data.center.weight = face->area();

  const size_t nfound = std::accumulate(region_found.begin(), region_found.end(), 0);
  if (nfound != region_found.size())
    throw std::runtime_error("fuck up my regions");

  return face_data;
}

}  // end namespace discretization

#endif
