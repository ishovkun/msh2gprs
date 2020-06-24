#include "IntegrationRule3dMS.hpp"
#include "../FeValues.hpp"  // provides FeValues

#ifdef WITH_EIGEN

namespace discretization {

IntegrationRule3dMS::IntegrationRule3dMS(PolyhedralElementBase & element)
    : _element(element)
{
  setup_storage_();

  const size_t n_parents = element._basis_functions.size();
  FeValues<angem::VTK_ID::TetrahedronID> fe_values;
  const auto & grid = element._element_grid;
  const Point parent_center = element._parent_cell.center();

  size_t icell = 0;
  for( auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell, ++icell  )
  {
    fe_values.update(*cell);
    const std::vector<size_t> & cell_verts = cell->vertices();
    const size_t nv = cell_verts.size();
    // auto & data = element._cell_data.points[cell->index()];
    auto & data = element._cell_data.points[icell];
    for (size_t parent_vertex=0; parent_vertex<n_parents; ++parent_vertex)
      for (size_t v=0; v<nv; ++v)
      {
        data.values[parent_vertex] += fe_values.value( v, 0 ) *
            element._basis_functions[parent_vertex][cell_verts[v]];
        data.grads[parent_vertex] += fe_values.grad(v, 0) *
            element._basis_functions[parent_vertex][cell_verts[v]];
      }
    data.weight = cell->volume();

    if (cell->polyhedron()->point_inside(parent_center))
    {
      _element._cell_data.center.values = data.values;
      _element._cell_data.center.grads = data.grads;
    }
  }
  _element._cell_data.center.weight = element._parent_cell.volume();
}

void IntegrationRule3dMS::setup_storage_()
{
  auto & data = _element._cell_data;
  const auto & basis_functions = _element._basis_functions;
  const size_t n_parent_vertices = _element._parent_cell.vertices().size();
  data.points.resize(_element._element_grid.n_active_cells());
  for (auto & point : data.points)
  {
    point.values.resize( basis_functions.size(), 0 );
    point.grads.resize( basis_functions.size() );
  }

  data.center.values.resize( basis_functions.size(), 0 );
  data.center.grads.resize( basis_functions.size() );
}

}  // end namespace discretization

#endif
