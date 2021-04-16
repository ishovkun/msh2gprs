#include "PolyhedralElementScaled.hpp"

namespace discretization {

using Point = angem::Point<3,double>;

PolyhedralElementScaled::PolyhedralElementScaled(const mesh::Cell & cell,
                                                 const mesh::Mesh & parent_grid,
                                                 PolyhedralElementBase & master,
                                                 const FiniteElementConfig & config)
    : PolyhedralElementBase(cell, parent_grid, config, true),
      _master(master)
{}

void PolyhedralElementScaled::build_fe_cell_data_()
{
  FiniteElementData const & master_data = _master.get_cell_data();
  size_t const nv = _parent_cell.n_vertices();
  size_t const nq = master_data.points.size();
  _cell_data.resize(nv, nq);

  auto const vert_coord = _parent_cell.vertex_coordinates();
  angem::Tensor2<3, double> du_dx;
  for (size_t q = 0; q < nq; ++q) {
    build_fe_point_data_(vert_coord, master_data.points[q],
                         _cell_data.points[q], du_dx);
  }

  build_fe_point_data_(vert_coord, master_data.center,
                       _cell_data.center, du_dx);
}

void PolyhedralElementScaled::
build_fe_point_data_(std::vector<angem::Point<3,double>> const & vertex_coord,
                     FEPointData const & master,
                     FEPointData & target,
                     angem::Tensor2<3, double> & du_dx) const
{
  double detJ;
  compute_detJ_and_invert_cell_jacobian_(master.grads,
                                         du_dx, detJ,
                                         _parent_cell.vertex_coordinates());

  target.weight = detJ * master.weight;
  target.values = master.values;

  if (detJ <= 0) throw std::runtime_error("Transformation det(J) is negative " + std::to_string(detJ));
  update_shape_grads_(master.grads, du_dx, target.grads);
}


void PolyhedralElementScaled::
update_shape_grads_(std::vector<angem::Point<3,double>> const & ref_grads,
                    angem::Tensor2<3, double> const & du_dx,
                    std::vector<angem::Point<3,double>> &grads) const
{
  for (size_t vertex = 0; vertex < grads.size(); ++vertex)
    grads[vertex] = {0.0, 0.0, 0.0};

  // compute the true shape function gradients
  for (size_t vertex = 0; vertex < grads.size(); ++vertex)
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        // d phi_vert / dx_i = (d phi_vert / d u_j) * (d u_j / d x_i)
        grads[vertex][i] += ref_grads[vertex][j] * du_dx(j, i);
}

void PolyhedralElementScaled::
compute_detJ_and_invert_cell_jacobian_(const std::vector<Point> & ref_grad,
                                       angem::Tensor2<3, double> & du_dx,
                                       double & detJ,
                                       std::vector<Point> const & vertex_coord) const
{
  // need to compute du/dx.
  // first compute dx/du = du/dξ * x
  angem::Tensor2<3, double> dx_du;
  for (size_t i=0; i<3; ++i)
    for (size_t j=0; j<3; ++j)
      for (size_t v=0; v < vertex_coord.size(); ++v)
        dx_du( i, j ) += ref_grad[v][j] * vertex_coord[v][i];

  // compute the determinant of transformation jacobian
  detJ = det(dx_du);
  // invert the jacobian to compute shape function gradients
  du_dx = invert(dx_du);
}

FiniteElementData PolyhedralElementScaled::
get_face_data(size_t iface, const angem::Basis<3,double> basis)
{
  auto const master_data = _master.get_face_data(iface, basis);
  size_t const nv = _parent_cell.n_vertices();
  size_t const nq = master_data.points.size();

  // check for the right numbering
  if ( nv != master_data.center.values.size())
    throw std::runtime_error("You need to renumber faces of the current polyhedron to match master");

  FiniteElementData face_data(nv, nq);
  auto const vert_coord = _parent_cell.faces()[iface]->vertex_coordinates();

  angem::Tensor2<3, double> du_dx;
  for (size_t q = 0; q < nq; ++q) {
    build_fe_point_data_(vert_coord, master_data.points[q],
                         face_data.points[q], du_dx);
  }

  build_fe_point_data_(vert_coord, master_data.center,
                       face_data.center, du_dx);

  return face_data;
}


}  // end namespace discretization
