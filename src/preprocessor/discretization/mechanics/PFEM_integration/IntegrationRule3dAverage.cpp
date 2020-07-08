#include "IntegrationRule3dAverage.hpp"
#include "../FeValues.hpp"

#ifdef WITH_EIGEN

namespace discretization {

using Point = angem::Point<3,double>;

IntegrationRule3dAverage::
IntegrationRule3dAverage(PolyhedralElementBase & element, const TributaryRegion3dBase  & tributary)
{
  setup_storage_(element, tributary);
  const auto & regions = tributary.get();
  const Point parent_center = element._parent_cell.center();
  bool center_found = false;
  // integrate over regions
  const size_t n_parents = element._basis_functions.size();
  std::vector<double> region_volumes( regions.size(), 0.0 );
  FeValues<angem::VTK_ID::TetrahedronID> fe_values;
  const auto & grid = element._subgrid;
  auto & cell_data = element._cell_data;
  for( auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell  )
  {
    const Point c = cell->center();
    fe_values.update(*cell);
    for (size_t region=0; region<regions.size(); ++region)  // tributary regions
    {
      if (regions[region].point_inside(c))
      {
        const std::vector<size_t> & cell_verts = cell->vertices();
        const size_t nv = cell_verts.size();
        region_volumes[region] += cell->volume();
        auto & data = cell_data.points[region];

        for (size_t parent_vertex=0; parent_vertex<n_parents; ++parent_vertex)
          for (size_t v=0; v<nv; ++v)
            for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
            {
              data.values[parent_vertex] += fe_values.value( v, q ) *
                  element._basis_functions[parent_vertex][cell_verts[v]] *
                  fe_values.JxW(q);
              data.grads[parent_vertex] += fe_values.grad(v, q) *
                  element._basis_functions[parent_vertex][cell_verts[v]] *
                  fe_values.JxW(q);
          }

        break;  // stop searching region
      }
    }
    const std::vector<size_t> & cell_verts = cell->vertices();
    const size_t nv = cell_verts.size();
    for (size_t parent_vertex=0; parent_vertex<n_parents; ++parent_vertex)
      for (size_t v=0; v<nv; ++v)
        for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
        {
          cell_data.center.values[parent_vertex] += fe_values.value(v, q) *
              element._basis_functions[parent_vertex][cell_verts[v]] *
              fe_values.JxW(q);
          cell_data.center.grads[parent_vertex] += fe_values.grad(v, q) *
              element._basis_functions[parent_vertex][cell_verts[v]] *
              fe_values.JxW(q);
        }
  }

  // normalize values and grads  by region volume
  for (size_t region=0; region<regions.size(); ++region)  // tributary regions
  {
    auto & data = cell_data.points[region];
    for (size_t parent_vertex=0; parent_vertex<n_parents; ++parent_vertex)
    {
      data.values[parent_vertex] /= region_volumes[region];
      data.grads[parent_vertex] /= region_volumes[region];
    }
    data.weight = region_volumes[region];
  }

  const double vol = element._parent_cell.volume();
  auto & data = cell_data.center;
  for (size_t parent_vertex=0; parent_vertex<n_parents; ++parent_vertex)
  {
    data.values[parent_vertex] /= vol;
    data.grads[parent_vertex] /= vol;
  }
  data.weight = vol;
}

void IntegrationRule3dAverage::setup_storage_(PolyhedralElementBase & element, const TributaryRegion3dBase  & tributary)
{
  auto & cell_data = element._cell_data;
  const auto & basis_functions = element._basis_functions;
  cell_data.points.resize( tributary.get().size() );
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
