#include "IntegrationRule3d.hpp"
#include "../FeValues.hpp"

namespace discretization {

IntegrationRule3d::IntegrationRule3d(PolyhedralElementBase const & element,
                                     TributaryRegion3dBase const  & tributary)
    : _nregions(tributary.size())
{
  size_t nq = 0;
  for (size_t region = 0; region < _nregions; ++region)
    nq += tributary.cells(region).size();
  nq += tributary.cells_center().size();
  _region.resize(nq);

  _data.resize(nq, FEPointData(element.n_vertices()));
  nq = 0;
  for (size_t region = 0; region < _nregions; ++region)
  {
    build_region_(element, tributary.cells(region), nq, region);
    nq += tributary.cells(region).size();
  }
  build_region_(element, tributary.cells_center(), nq, _nregions);
}

void IntegrationRule3d::build_region_(PolyhedralElementBase const & element,
                                      std::vector<size_t> const & cells,
                                      size_t offset,
                                      size_t region)
{
  size_t const npv = element.n_vertices();
  const size_t nv = ElementTraits<angem::TetrahedronID>::n_vertices;  // number of vertices in triangle
  FeValues<angem::VTK_ID::TetrahedronID> fe_values;
  const auto & grid = element.get_grid();
  auto const & sf = element.get_basis_functions();

  for (size_t icell = 0; icell < cells.size(); ++icell)
  {
    auto const & cell = grid.cell(cells[icell]);
    size_t qg = offset + icell;
    _region[qg] = region;
    const std::vector<size_t> & cell_verts = cell.vertices();
    fe_values.update(cell);
    for (size_t pv = 0; pv < npv; ++pv) {
      for (size_t v = 0; v < nv; ++v)
      {
        double const parent_shape_value = sf[pv][cell_verts[v]];
        for (size_t q = 0; q < fe_values.n_integration_points(); ++q) {
          _data[qg].values[pv] += fe_values.value(v, q) *
                                     parent_shape_value *
                                     fe_values.JxW(q);
          _data[qg].grads[pv] += fe_values.grad(v, q) *
                                    parent_shape_value *
                                    fe_values.JxW(q);
        }
      }
    }

    for (size_t q = 0; q < fe_values.n_integration_points(); ++q) {
      _data[qg].weight += fe_values.JxW(q);
    }
  }
}

FiniteElementData IntegrationRule3d::
integrate(std::vector<angem::Point<3,double>> const & vertices) const
{
  angem::Tensor2<3, double> du_dx;
  FiniteElementData result(vertices.size(), _nregions);
  for (size_t q = 0; q < _data.size(); ++q)
  {
    if (_region[q] < _nregions)
      build_fe_point_data_(vertices, _data[q], result.points[_region[q]], du_dx);
    else
      build_fe_point_data_(vertices, _data[q], result.center, du_dx);
  }
  return result;
}


void IntegrationRule3d::
build_fe_point_data_(std::vector<angem::Point<3,double>> const & vertex_coord,
                     FEPointData const & master,
                     FEPointData & target,
                     angem::Tensor2<3, double> & du_dx) const
{
  double const detJ = compute_detJ_and_invert_cell_jacobian_(master.grads, du_dx, vertex_coord);
  if (detJ <= 0)
    throw std::runtime_error("Cell Transformation det(J) is negative " + std::to_string(detJ));
  target.weight += detJ * master.weight;
  for (size_t v = 0; v < vertex_coord.size(); ++v)
  {
    target.values[v] += master.values[v] * detJ * master.weight;
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        target.grads[v][i] += master.grads[v][j] * du_dx(j, i) * detJ * master.weight;
  }
}

double IntegrationRule3d::
compute_detJ_and_invert_cell_jacobian_(const std::vector<Point> & ref_grad,
                                       angem::Tensor2<3, double> & du_dx,
                                       std::vector<Point> const & vertex_coord) const
{
  // compute the true shape function gradients
  // i, j - component indices
  // k - vertex index
  // ξ - coordinate in master element
  // x - coordinate in current basis
  // Ψ = Ψ(ξ) master element shape function

  // first compute ∂xᵢ/dξⱼ = Σⱼ∂Ψₖ/∂ξⱼ * xₖᵢ
  // xₖᵢ - i-component of kth vertex coordinate
  angem::Tensor2<3, double> dx_du;
  for (size_t i=0; i<3; ++i)
    for (size_t j=0; j<3; ++j)
      for (size_t v=0; v < vertex_coord.size(); ++v)
        dx_du( i, j ) += ref_grad[v][j] * vertex_coord[v][i];

  // compute the determinant of transformation jacobian
  double const detJ = det(dx_du);
  // invert the jacobian to compute shape function gradients ∂ξᵢ/∂xⱼ = (∂xᵢ/dξⱼ)⁻¹
  du_dx = invert(dx_du);
  return detJ;
}


}  // end namespace discretization
