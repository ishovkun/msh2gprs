#include "IntegrationRule2d.hpp"
#include "../FeValues.hpp"

namespace discretization {

IntegrationRule2d::IntegrationRule2d(const PolyhedralElementBase & element,
                                     const TributaryRegion2dBase  & tributary,
                                     const size_t parent_face,
                                     const angem::Basis<3, double> & basis)
    : _nregions(tributary.size())
{
  compute_parent_vertices_(element.host_cell(), parent_face);

  size_t nq = 0;
  for (size_t region = 0; region < _nregions; ++region)
    nq += tributary.faces(region).size();
  nq += tributary.faces_center().size();
  _region.resize(nq);

  _data.resize(nq, FEPointData(_npv));
  // then fill the storage
  nq = 0;
  for (size_t region = 0; region < _nregions; ++region)
  {
    build_region_(element, tributary.faces(region), basis, nq, region);
    nq += tributary.faces(region).size();
  }
  build_region_(element, tributary.faces_center(), basis, nq, _nregions);
}

void IntegrationRule2d::build_region_(PolyhedralElementBase const & element,
                                      std::vector<size_t> const & faces,
                                      const angem::Basis<3, double> & basis,
                                      size_t offset, size_t region)
{
  // Loop over all triangular faces of the subgrid and save they weight, values, and grads
  auto const & grid = element.get_grid();
  FeValues<angem::VTK_ID::TriangleID> fe_values;
  fe_values.set_basis(basis);
  auto const & sf = element.get_basis_functions();
  const size_t nv = ElementTraits<angem::TriangleID>::n_vertices;  // number of vertices in triangle

  for (size_t iface = 0; iface < faces.size(); ++iface)
  {
    const mesh::Face & face = grid.face(faces[iface]);
    size_t const qg = offset + iface;
    _region[qg] = region;
    fe_values.update(face);
    const std::vector<size_t> & face_verts = face.vertices();

    for (size_t q = 0; q < fe_values.n_integration_points(); ++q)
      _data[qg].weight += fe_values.JxW(q);

    for (size_t pv = 0; pv < _npv; ++pv) {
      for (size_t v = 0; v < nv; ++v) {
        double const parent_shape_value = sf[_parent_vertices[pv]][face_verts[v]];
        for (size_t q = 0; q < fe_values.n_integration_points(); ++q) {
          _data[qg].values[pv] += fe_values.value(v, q) * parent_shape_value;
          _data[qg].grads[pv] += fe_values.grad(v, q) * parent_shape_value;
        }
      }
    }
  }
}

FiniteElementData IntegrationRule2d::
integrate(std::vector<angem::Point<3,double>> const & vertices,
          angem::Basis<3, double>             const & basis) const
{
  assert( _npv == vertices.size() && "number of vertices does not match" );

  angem::Tensor2<3, double> du_dx;
  FiniteElementData result(vertices.size(), _nregions);
  for (size_t q = 0; q < _data.size(); ++q) {
    if (_region[q] < _nregions)
      build_fe_point_data_(vertices, _data[q], result.points[_region[q]], basis, du_dx);
    else
      build_fe_point_data_(vertices, _data[q], result.center,             basis, du_dx);
  }

  for (size_t region = 0; region < _nregions; ++region)
    scale_data_(result.points[region]);
  scale_data_(result.center);

  return result;
}

void IntegrationRule2d::build_fe_point_data_(std::vector<angem::Point<3,double>> const & vertex_coord,
                                             FEPointData const & master,
                                             FEPointData & target,
                                             const angem::Basis<3, double> & basis,
                                             angem::Tensor2<3, double> & du_dx) const
{
  double const detJ = compute_inverse_jacobian_(master.grads, vertex_coord, basis, du_dx);
  if ( detJ <= 0 ) throw std::runtime_error("Face Transformation det(J) is negative " + std::to_string(detJ));

  target.weight += detJ * master.weight;
  for (size_t v = 0; v < _npv; ++v)
  {
    target.values[v] += master.values[v] * detJ * master.weight;
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        target.grads[v][i] += master.grads[v][j] * du_dx(j, i) * detJ * master.weight;
  }

}

double IntegrationRule2d::compute_inverse_jacobian_(std::vector<angem::Point<3,double>> const & ref_grad,
                                                    std::vector<angem::Point<3,double>> const & vertex_coord,
                                                    angem::Basis<3, double>             const & basis,
                                                    angem::Tensor2<3, double>                 & J_inv) const
{
  angem::Plane<double> plane(vertex_coord);
  plane.set_basis(basis);
  std::vector<angem::Point<3,double>> loc_coord(_npv);
  for (size_t v = 0; v < _npv; ++v)
  {
    loc_coord[v] = plane.local_coordinates(vertex_coord[v]);
    // std::cout << "loc coord (" << v << ") = " << loc_coord[v] << std::endl;
    if ( std::fabs(loc_coord[v][2]) > 1e-10 )
      logging::warning() << "non-planar face or basis is set for non-planar surfaces" << std::endl;
  }

  angem::Tensor2<3, double> J;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      for (size_t v = 0; v < _npv; ++v)
        J(i, j) += ref_grad[v][j] * loc_coord[v][i];
  J(2, 2) = 1;
  J_inv = invert(J);
  // std::cout << "ref grad: ";
  // for (size_t v = 0; v <  _npv; ++v)
  //   std::cout << v << " " << ref_grad[v] << std::endl;

  // std::cout << "J = " << J << std::endl;
  // std::cout << "J_1 = " << J_inv << std::endl;
  return det(J);
}

void IntegrationRule2d::compute_parent_vertices_(mesh::Cell const & cell, size_t iface)
{
  auto faces = cell.faces();
  _npv = faces[iface]->vertices().size();
  const auto cv = cell.vertices();
  const auto fv = faces[iface]->vertices();
  std::unordered_map<size_t, size_t> mp;
  for (size_t v = 0; v < cv.size(); ++v)
    mp[cv[v]] = v;

  _parent_vertices.resize(_npv);
  for (size_t v = 0; v < _npv; ++v)
    _parent_vertices[v] = mp[ fv[v] ];

  // const auto face = _element._parent_cell.faces()[_parent_face];
    // for (size_t v=0; v < npv; ++v)
    //   _parent_vertices[v] =
    //       std::distance(parent_cell_vertices.begin(),
    //                     std::find( parent_cell_vertices.begin(), parent_cell_vertices.end(),
    //                                face->vertices()[v]));
}

void IntegrationRule2d::scale_data_(FEPointData & data) const
{
  for (size_t pv = 0; pv < _npv; ++pv)
  {
    data.values[pv] /= data.weight;
    data.grads[pv] /= data.weight;
  }
}

}  // end namespace discretization
