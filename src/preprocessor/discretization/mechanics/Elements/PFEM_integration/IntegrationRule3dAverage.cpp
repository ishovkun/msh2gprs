#include "IntegrationRule3dAverage.hpp"
#include "../FeValues.hpp"

#ifdef WITH_EIGEN

namespace discretization {

using Point = angem::Point<3,double>;

IntegrationRule3dAverage::
IntegrationRule3dAverage(PolyhedralElementBase & element, const TributaryRegion3dBase  & tributary)
    : _element(element)
{
  // setup storage
  element._cell_data.resize( element._basis_functions.size(), tributary.size() );
  // compute fe data in gauss points
  for (size_t iregion = 0; iregion < tributary.size(); ++iregion)
    compute_fe_values_(tributary.cells(iregion), _element._cell_data.points[iregion]);
  // compute fe data center
  compute_fe_values_(tributary.cells_center(), _element._cell_data.center);
}

void IntegrationRule3dAverage::compute_fe_values_(const std::vector<size_t> &cells,
                                                  FEPointData & dst)
{
  const size_t n_parent_vertices = _element._basis_functions.size();
  const size_t nv = ElementTraits<angem::TetrahedronID>::n_vertices;  // number of vertices in triangle
  FeValues<angem::VTK_ID::TetrahedronID> fe_values;
  const auto & grid = _element._subgrid;
  double sum_volumes = 0.f;
  for (const size_t icell : cells) {
    const auto & cell = grid.cell(icell);
    const std::vector<size_t> & cell_verts = cell.vertices();
    fe_values.update(cell);

    for (size_t parent_vertex = 0; parent_vertex < n_parent_vertices; ++parent_vertex) {
      for (size_t v = 0; v < nv; ++v) {
        double const parent_shape_value = _element._basis_functions[parent_vertex][cell_verts[v]];
        for (size_t q = 0; q < fe_values.n_integration_points(); ++q) {
          dst.values[parent_vertex] += fe_values.value(v, q) * parent_shape_value * fe_values.JxW(q);
          dst.grads[parent_vertex]  += fe_values.grad(v, q)  * parent_shape_value * fe_values.JxW(q);
        }
      }
    }

    for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
      sum_volumes += fe_values.JxW(q);
  }

  for (size_t parent_vertex=0; parent_vertex < n_parent_vertices; ++parent_vertex) {
    dst.values[parent_vertex] /= sum_volumes;
    dst.grads[parent_vertex] /= sum_volumes;
  }

  dst.weight = sum_volumes;

  if (sum_volumes < 1e-10)
    throw std::runtime_error("Messed up volumes");
}

}  // end namespace discretization


#endif
