#include "IntegrationRule3dAverage.hpp"
#include "../FeValues.hpp"

#ifdef WITH_EIGEN

namespace discretization {

using Point = angem::Point<3,double>;

IntegrationRule3dAverage::
IntegrationRule3dAverage(PolyhedralElementBase & element, const TributaryRegion3dBase  & tributary)
    : _element(element)
{
  setup_storage_(element, tributary);
  for (size_t iregion = 0; iregion < tributary.size(); ++iregion)
  {
    compute_fe_values_(tributary.get(iregion), _element._cell_data.points[iregion]);
    _element._cell_data.points[iregion].weight = tributary.volume(iregion);
  }

  compute_fe_values_(tributary.get_center(), _element._cell_data.center);
  _element._cell_data.center.weight = tributary.volume_center();
}

void IntegrationRule3dAverage::compute_fe_values_(const std::vector<size_t> &cells,
                                                  FEPointData & dst)
{
  const size_t n_parents = _element._basis_functions.size();
  FeValues<angem::VTK_ID::TetrahedronID> fe_values;
  const auto & grid = _element._subgrid;
  double sum_volumes = 0;
  for (const size_t icell : cells)
  {
    const auto & cell = grid.cell(icell);
    const std::vector<size_t> & cell_verts = cell.vertices();
    const size_t nv = cell_verts.size();
    fe_values.update(cell);
    for (size_t parent_vertex=0; parent_vertex<n_parents; ++parent_vertex)
      for (size_t v=0; v<nv; ++v)
        for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
        {
          dst.values[parent_vertex] += fe_values.value( v, q ) *
              _element._basis_functions[parent_vertex][cell_verts[v]] *
              fe_values.JxW(q);
          dst.grads[parent_vertex] += fe_values.grad(v, q) *
              _element._basis_functions[parent_vertex][cell_verts[v]] *
              fe_values.JxW(q);
        }

    for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
      sum_volumes += fe_values.JxW(q);
  }

  for (size_t parent_vertex=0; parent_vertex<n_parents; ++parent_vertex)
  {
    dst.values[parent_vertex] /= sum_volumes;
    dst.grads[parent_vertex] /= sum_volumes;
  }

  if (sum_volumes < 1e-10)
    throw std::runtime_error("Messed up volumes");
}

void IntegrationRule3dAverage::setup_storage_(PolyhedralElementBase & element, const TributaryRegion3dBase  & tributary)
{
  auto & cell_data = element._cell_data;
  const auto & basis_functions = element._basis_functions;
  cell_data.points.resize( tributary.size() );
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
