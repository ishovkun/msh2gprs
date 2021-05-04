#include "PolyhedralElementScaled.hpp"
#include "logger/Logger.hpp"
// debug
// #include "PFEM_integration/TributaryRegion3dFaces.hpp"
// #include "PFEM_integration/IntegrationRule3dAverage.hpp"
#include "PFEM_integration/IntegrationRule3d.hpp"

namespace discretization {

using Point = angem::Point<3,double>;

PolyhedralElementScaled::PolyhedralElementScaled(const mesh::Cell & cell,
                                                 const mesh::Mesh & parent_grid,
                                                 PolyhedralElementBase & master,
                                                 const FiniteElementConfig & config)
    : PolyhedralElementBase(cell, parent_grid, config, true),
      _master(master)
{
  map_vertices_to_master_();
}

void PolyhedralElementScaled::map_vertices_to_master_()
{
  auto const & v_cur = _parent_cell.vertices();
  auto const & v_master = _master.host_cell().vertices();
  for (size_t v = 0; v < _parent_cell.n_vertices(); ++v)
    _vertex_mapping[v_cur[v]] = v_master[v];
}

void PolyhedralElementScaled::build_fe_cell_data_()
{
  // size_t const npv = _parent_cell.n_vertices();
  auto const vert_coord = _parent_cell.vertex_coordinates();
  auto const & rule = _master.integration_rule3();
  _cell_data = rule.integrate(_parent_cell.vertex_coordinates());

  // // new shit
  // _subgrid =  _master.get_grid();
  // _basis_functions = _master.get_basis_functions();

  // // move vertices
  // for (size_t v = 0; v < _subgrid.n_vertices(); ++v)
  // {
  //    _subgrid.vertex(v).set_zero();
  //   for (size_t pv = 0; pv < npv; ++pv) {
  //     _subgrid.vertex(v) += vert_coord[pv] * _basis_functions[pv][v];
  //   }
  // }

  // TributaryRegion3dFaces tributary3d(*this);
  // IntegrationRule3dAverage rule_cell(*this, tributary3d);

  // version 2
  // IntegrationRule3dFull rule(_master);
  // TributaryRegion3dFaces tributary3d(_master);
  // size_t nq = tributary3d.size();
  // FiniteElementData const & master_data = _master.get_cell_data();
  // _cell_data.resize(npv, nq);
  // angem::Tensor2<3, double> du_dx;
  // for (size_t q = 0; q < nq; ++q)
  // {
  //   // std::cout << "q = " << q << std::endl;
  //   // std::cout << "master_data.points.size() = " << master_data.points.size() << std::endl;
  //   // std::cout << "_cell_data.points.size() = " << _cell_data.points.size() << std::endl;
  //   // std::cout << "indices: ";
  //   // for (size_t icell : tributary3d.get_indices(q))
  //   //   std::cout << icell << " ";
  //   // std::cout << std::endl;
  //   for (size_t icell : tributary3d.cells(q))
  //   {
  //     // std::cout << "icell = " << icell << std::endl;
  //     icell = icell - 1;  // active cell indices shifted by 1
  //     build_fe_point_data_append_(vert_coord, master_data.points[icell], _cell_data.points[q], du_dx);
  //     build_fe_point_data_append_(vert_coord, master_data.points[icell], _cell_data.center, du_dx);
  //   }

  //   // scale
  //   for (size_t pv = 0; pv < npv; ++pv)
  //   {
  //    _cell_data.points[q].values[pv] /= _cell_data.points[q].weight;
  //    _cell_data.points[q].grads[pv] /= _cell_data.points[q].weight;
  //   }
  // }

    // scale
    // for (size_t pv = 0; pv < npv; ++pv)
    // {
    //  _cell_data.center.values[pv] /= _cell_data.center.weight;
    //  _cell_data.center.grads[pv] /= _cell_data.center.weight;
    // }

  // angem::Tensor2<3, double> du_dx;
  // for (size_t q = 0; q < nq; ++q) {
  //   build_fe_point_data_(vert_coord, master_data.points[q], _cell_data.points[q], du_dx);
  // }

  // build_fe_point_data_(vert_coord, master_data.center, _cell_data.center, du_dx);
}

void PolyhedralElementScaled::
build_fe_point_data_append_(std::vector<angem::Point<3,double>> const & vertex_coord,
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

void PolyhedralElementScaled::
build_fe_point_data_(std::vector<angem::Point<3,double>> const & vertex_coord,
                     FEPointData const & master,
                     FEPointData & target,
                     angem::Tensor2<3, double> & du_dx) const
{
  double const detJ = compute_detJ_and_invert_cell_jacobian_(master.grads, du_dx, vertex_coord);
  if (detJ <= 0) throw std::runtime_error("Cell Transformation det(J) is negative " + std::to_string(detJ));

  target.weight = detJ * master.weight;
  target.values = master.values;
  update_shape_grads_(master.grads, du_dx, target.grads);
}


void PolyhedralElementScaled::
update_shape_grads_(std::vector<angem::Point<3,double>> const & ref_grads,
                    angem::Tensor2<3, double> const & du_dx,
                    std::vector<angem::Point<3,double>> &grads) const
{
  // compute the true shape function gradients
  // i, j - component indices
  // k - vertex index
  // ξ - coordinate in master element
  // x - coordinate in current basis
  // φ = φ(x) current element shape function
  // Ψ = Ψ(ξ) master element shape function
  // ∂φₖ / ∂xᵢ = Σⱼ (∂Ψₖ / d ξⱼ) * (∂ξⱼ / dxᵢ)
  for (size_t vertex = 0; vertex < grads.size(); ++vertex)
    grads[vertex].set_zero();
  for (size_t k = 0; k < grads.size(); ++k)
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        grads[k][i] += ref_grads[k][j] * du_dx(j, i);
}

double PolyhedralElementScaled::
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

FiniteElementData PolyhedralElementScaled::get_face_data(size_t iface)
{
  auto const * const face = _parent_cell.faces()[iface];
  auto const * const master_face = _master.host_cell().faces()[iface];
  auto vert_coord = _parent_cell.faces()[iface]->vertex_coordinates();
  reorder_face_vertices_(*face, *master_face, vert_coord);

  auto basis = get_face_basis_(*face, host_cell());
  FiniteElementData const master_data = _master.get_face_data(iface);
  size_t const nv = vert_coord.size();
  size_t const nq = master_data.points.size();

  // check for the right numbering
  if ( nv != master_data.center.values.size() )
    throw std::runtime_error("You need to renumber faces of the current polyhedron to match master");

  FiniteElementData face_data(nv, nq);

   angem::Tensor2<3, double> du_dx;
   for (size_t q = 0; q < nq; ++q)
     build_fe_face_point_data_(vert_coord, master_data.points[q],
                               face_data.points[q], du_dx, basis);

   build_fe_face_point_data_(vert_coord, master_data.center, face_data.center,
                             du_dx, basis);
   // FIXME: we reordered vertices, so should probably reorder the resulting values back
   // if ( reverse )
   // {
   //   for (size_t q = 0; q < nq; ++q)
   //   {
   //     std::reverse( face_data.points[q].values.begin(),
   //                   face_data.points[q].values.end());
   //     std::reverse( face_data.points[q].grads.begin(),
   //                   face_data.points[q].grads.end());
   //   }
   //   std::reverse( face_data.center.values.begin(),
   //                 face_data.center.values.end());
   //   std::reverse( face_data.center.grads.begin(),
   //                 face_data.center.grads.end());
   // }
   return face_data;
}

void PolyhedralElementScaled::
build_fe_face_point_data_(std::vector<angem::Point<3,double>> const & vertex_coord,
                          FEPointData const & master,
                          FEPointData & current,
                          angem::Tensor2<3, double> & du_dx,
                          const angem::Basis<3,double> & basis) const
{

  double const detJ = compute_detJ_and_invert_face_jacobian_(master.grads, du_dx, vertex_coord, basis);
  if ( detJ <= 0 ) throw std::runtime_error("Face Transformation det(J) is negative " + std::to_string(detJ));
  current.values = master.values;
  current.weight = master.weight * detJ;
  update_shape_grads_(master.grads, du_dx, current.grads);
}

double PolyhedralElementScaled::
compute_detJ_and_invert_face_jacobian_(const std::vector<angem::Point<3,double>> & ref_grad,
                                       angem::Tensor2<3, double> & J_inv,
                                       std::vector<angem::Point<3,double>> const & vertex_coord,
                                       const angem::Basis<3,double> & basis) const
{
  angem::Plane<double> plane(vertex_coord);
  plane.set_basis(basis);
  size_t const nv = vertex_coord.size();
  std::vector<Point> loc_coord(nv);
  for (size_t v = 0; v < nv; ++v)
  {
    loc_coord[v] = plane.local_coordinates(vertex_coord[v]);
    if ( std::fabs(loc_coord[v][2]) > 1e-10 )
      logging::warning() << "non-planar face or basis is set for non-planar surfaces" << std::endl;
  }

  angem::Tensor2<3, double> J;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      for (size_t v = 0; v < nv; ++v)
        J(i, j) += ref_grad[v][j] * loc_coord[v][i];
  J(2, 2) = 1;
  J_inv = invert(J);
  return det(J);
}

void PolyhedralElementScaled::
reorder_face_vertices_(mesh::Face const & target,
                       mesh::Face const & master,
                       std::vector<angem::Point<3,double>> & coords)
{
  auto const & verts = target.vertices();
  auto const & master_verts = master.vertices();
  std::vector<size_t> order(verts.size(), 0u);
  for (size_t v = 0; v < verts.size(); ++v)
  {
    size_t const mapped = _vertex_mapping[verts[v]];
    auto it = std::find( master_verts.begin(), master_verts.end(), mapped);
    if (it == master_verts.end())  // just debugging
    {
      std::cout << "cannot find mapped vertex " << verts[v] << " (mapped as "
                << _vertex_mapping[ verts[v] ] << ")" << std::endl;
      std::cout << "mapping: " << std::endl;
      for (auto it1 : _vertex_mapping)
          std::cout << it1.first << "-" << it1.second << std::endl;
      std::cout << "master verts" << std::endl;
      for (auto v : master_verts)
        std::cout << v << " ";
      std::cout << std::endl;
      throw std::runtime_error("something's wrong I can feel it");
    }

    order[v] = std::distance(master_verts.begin(), it);
  }

  auto const copy = coords;
  for (size_t v = 0; v < verts.size(); ++v)
    coords[v] = copy[order[v]];
}

angem::Polyhedron<double>
PolyhedralElementScaled::create_pyramid_(const std::vector<size_t> & face,
                                         const std::vector<Point> & vertices) const
{
  const size_t vertex_center = vertices.size() - 1;  // HACK: I just pushed it to this array
  const auto c = _parent_cell.center();
  std::vector<std::vector<size_t>> pyramid_faces;
  for (size_t iv=0; iv<face.size(); ++iv)
  {
    size_t const v1 = face[iv];
    size_t const v2 = face[(iv+1) % face.size()];
    pyramid_faces.push_back( {v1, v2, vertex_center} );
  }
  pyramid_faces.push_back( face );  // base

  angem::Polyhedron<double> pyramid(vertices, pyramid_faces);
  return pyramid;
}


}  // end namespace discretization
