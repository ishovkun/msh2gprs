#include "IntegrationRule3dPointwise.hpp"
#include "../FeValues.hpp"

#ifdef WITH_EIGEN

namespace discretization {

using Point = angem::Point<3,double>;

IntegrationRule3dPointwise::IntegrationRule3dPointwise(PolyhedralElementBase & element, const TributaryRegion3dBase  & tributary)
    : _element(element), _tributary(tributary)
{
  setup_storage_();
  const auto & regions = tributary.get();
  const Point parent_center = element._parent_cell.center();
  bool center_found = false;
  // integrate over regions
  const size_t n_parents = element._basis_functions.size();
  FeValues<angem::VTK_ID::TetrahedronID> fe_values;
  const auto & grid = element._element_grid;
  auto & cell_data = element._cell_data;
  for( auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell  )
  {
    for (size_t region=0; region<regions.size(); ++region)  // tributary regions
    {
      const Point c = regions[region].center();
      // if (regions[region].point_inside(c))
      if (cell->polyhedron()->point_inside(c))
      {
        fe_values.update(*cell, {c});
        const std::vector<size_t> & cell_verts = cell->vertices();
        const size_t nv = cell_verts.size();
        auto & data = cell_data.points[region];
        const size_t q = 0;

        for (size_t parent_vertex=0; parent_vertex<n_parents; ++parent_vertex)
          for (size_t v=0; v<nv; ++v)
            {
              data.values[parent_vertex] += fe_values.value( v, q ) *
                  element._basis_functions[parent_vertex][cell_verts[v]];
              data.grads[parent_vertex] += fe_values.grad(v, q) *
                  element._basis_functions[parent_vertex][cell_verts[v]];
          }

        break;  // stop searching region
      }
    }

    if (!center_found)
      if (cell->polyhedron()->point_inside(parent_center))
      {
        center_found = true;
        fe_values.update(*cell, {parent_center});

        const std::vector<size_t> & cell_verts = cell->vertices();
        const size_t nv = cell_verts.size();
        auto & data = cell_data.center;
        const size_t q = 0;

        for (size_t parent_vertex=0; parent_vertex<n_parents; ++parent_vertex)
          for (size_t v=0; v<nv; ++v)
            {
              data.values[parent_vertex] += fe_values.value( v, q ) *
                  element._basis_functions[parent_vertex][cell_verts[v]];
              data.grads[parent_vertex] += fe_values.grad(v, q) *
                  element._basis_functions[parent_vertex][cell_verts[v]];
          }
        // throw std::runtime_error("not implemented");
       
      }
  }

  // normalize values and grads  by region volume
  for (size_t region=0; region<regions.size(); ++region)  // tributary regions
  {
    auto & data = cell_data.points[region];
    data.weight = regions[region].volume();
  }
}

void IntegrationRule3dPointwise::setup_storage_()
{
  auto & cell_data = _element._cell_data;
  const auto & basis_functions = _element._basis_functions;
  cell_data.points.resize( _tributary.get().size() );
  for (auto & point : cell_data.points)
  {
    point.values.resize( basis_functions.size(), 0 );
    point.grads.resize( basis_functions.size() );
  }

  cell_data.center.values.resize( basis_functions.size(), 0 );
  cell_data.center.grads.resize( basis_functions.size() );
}

}  // end namespace discretization

#endif
