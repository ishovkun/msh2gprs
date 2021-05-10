#include "IntegrationRuleFracture.hpp"
#include "../FeValues.hpp"

namespace discretization {

IntegrationRuleFracture::
IntegrationRuleFracture(PolyhedralElementBase const & element,
                        TributaryRegion2dBase const & tributary,
                        size_t                        parent_face,
                        angem::Basis<3, double> const & basis)
    : _nregions(tributary.size())
{
  compute_parent_vertices_(element.host_cell(), parent_face);

  // setup storage
  size_t nq = 0;
  for (size_t region = 0; region < _nregions; ++region)
    nq += tributary.faces(region).size();
  nq += tributary.faces_center().size();
  _region.resize(nq);
  _data.resize(nq, FEPointData(_npcv));

  // fill storage
  nq = 0;  // offset
  for (size_t region = 0; region < _nregions; ++region)
  {
    build_region_(element, tributary.faces(region), basis, nq, region);
    nq += tributary.faces(region).size();
  }
  build_region_(element, tributary.faces_center(), basis, nq, _nregions);
}

void IntegrationRuleFracture::compute_parent_vertices_(mesh::Cell const & cell, size_t iface)
{
  _npcv = cell.n_vertices();
  auto faces = cell.faces();
  _npfv = faces[iface]->vertices().size();
  const auto cv = cell.vertices();
  const auto fv = faces[iface]->vertices();
  std::unordered_map<size_t, size_t> mp;
  for (size_t v = 0; v < cv.size(); ++v)
    mp[cv[v]] = v;

  _parent_vertices.resize(_npfv);
  for (size_t v = 0; v < _npfv; ++v)
    _parent_vertices[v] = mp[ fv[v] ];
}

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

void IntegrationRuleFracture::build_region_(PolyhedralElementBase const & element,
                                            std::vector<size_t> const & faces,
                                            const angem::Basis<3, double> & basis,
                                            size_t offset, size_t region)
{
  /* Special case of face integration rule:
   * Need values of all cell shape functions within face integration
   * points. */
  auto const & grid = element.get_grid();
  FeValues<angem::VTK_ID::TetrahedronID> fe_cell_values;
  FeValues<angem::VTK_ID::TriangleID> fe_face_values;
  fe_face_values.set_basis(basis);
  auto const & sf = element.get_basis_functions();
  const size_t nv = ElementTraits<angem::TetrahedronID>::n_vertices;
  std::vector<Point> local_integration_points;

  for (size_t iface = 0; iface < faces.size(); ++iface)
  {
    mesh::Face const & face = grid.face(faces[iface]);
    mesh::Cell const & cell = *face.neighbors()[0];  // face only has one neighbor
    std::vector<size_t> const & cell_verts = cell.vertices();
    size_t const qg = offset + iface;

    fe_face_values.update(face);
    get_face_integration_points(fe_face_values, local_integration_points);
    fe_cell_values.update(cell, local_integration_points);

    for (size_t q = 0; q < fe_face_values.n_integration_points(); ++q)
      _data[qg].weight += fe_face_values.JxW(q);

    for (size_t pv = 0; pv < _npcv; ++pv) {
      for (size_t v = 0; v < nv; ++v) {
        double const parent_shape_value = sf[pv][cell_verts[v]];
        for (size_t q = 0; q < fe_cell_values.n_integration_points(); ++q) {
          _data[qg].values[pv] += fe_cell_values.value(v, q) * parent_shape_value;
          _data[qg].grads[pv]  += fe_cell_values.grad(v, q)  * parent_shape_value;
        }
      }
    }

  }
}

FiniteElementData IntegrationRuleFracture::integrate(std::vector<angem::Point<3,double>> const & vertices,
                                                     angem::Basis<3, double>             const & basis) const
{
  assert( _npcv == vertices.size() && "number of vertices does not match" );

  FiniteElementData result(vertices.size(), _nregions);
  for (size_t q = 0; q < _data.size(); ++q) {
    if (_region[q] < _nregions)
      build_fe_point_data_(vertices, _data[q], result.points[_region[q]], basis);
    else
      build_fe_point_data_(vertices, _data[q], result.center,             basis);
  }

  for (size_t region = 0; region < _nregions; ++region)
    scale_data_(result.points[region]);
  scale_data_(result.center);

  return result;
}

void IntegrationRuleFracture::build_fe_point_data_(std::vector<angem::Point<3,double>> const & vertex_coord,
                                                   FEPointData const & master,
                                                   FEPointData & target,
                                                   const angem::Basis<3, double> & basis) const
{
  // compute face transformation jacobian to get the weights
  double const detJ_face = compute_face_scaling_(master.grads, vertex_coord, basis);

  // compute cell transformation to compute values and grads
  angem::Tensor2<3, double> du_dx;
  double const detJ = compute_inverse_cell_jacobian_(master.grads, du_dx, vertex_coord);
  if ( detJ <= 0 ) throw std::runtime_error("Face Transformation det(J) is negative " + std::to_string(detJ));

  target.weight += detJ_face * master.weight;
  for (size_t v = 0; v < _npcv; ++v)
  {
    target.values[v] += master.values[v] * detJ_face * master.weight;
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        target.grads[v][i] += master.grads[v][j] * du_dx(j, i) * detJ_face * master.weight;
  }
}

double IntegrationRuleFracture::compute_face_scaling_(std::vector<angem::Point<3,double>> const & ref_grad,
                                                      std::vector<angem::Point<3,double>> const & vertex_coord,
                                                      angem::Basis<3, double>             const & basis) const
{
  angem::Plane<double> plane(vertex_coord[_parent_vertices[0]], basis(2));
  plane.set_basis(basis);
  std::vector<angem::Point<3,double>> loc_coord(_npfv);
  for (size_t v = 0; v < _npfv; ++v)
  {
    loc_coord[v] = plane.local_coordinates(vertex_coord[_parent_vertices[v]]);
    if ( std::fabs(loc_coord[v][2]) > 1e-10 )
      logging::warning() << "non-planar face or basis is set for non-planar surfaces" << std::endl;
  }

  angem::Tensor2<3, double> J;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      for (size_t v = 0; v < _npfv; ++v)
        J(i, j) += ref_grad[v][j] * loc_coord[v][i];
  J(2, 2) = 1;
  return det(J);
}

double IntegrationRuleFracture::
compute_inverse_cell_jacobian_(std::vector<angem::Point<3,double>> const & ref_grad,
                               angem::Tensor2<3, double> & du_dx,
                               std::vector<angem::Point<3,double>> const & vertex_coord) const
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
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t v = 0; v < vertex_coord.size(); ++v)
        dx_du( i, j ) += ref_grad[v][j] * vertex_coord[v][i];

  // compute the determinant of transformation jacobian
  double const detJ = det(dx_du);
  // invert the jacobian to compute shape function gradients ∂ξᵢ/∂xⱼ = (∂xᵢ/dξⱼ)⁻¹
  du_dx = invert(dx_du);
  return detJ;
}

void IntegrationRuleFracture::scale_data_(FEPointData & data) const
{
  for (size_t pv = 0; pv < _npcv; ++pv)
  {
    data.values[pv] /= data.weight;
    data.grads[pv]  /= data.weight;
  }
}


}  // end namespace discretization
