#include "DiscretizationPolyhedralFEMOptimized.hpp"
#include "Isomorphism.hpp"         // provides isomorphism


namespace discretization {
using Point = angem::Point<3,double>;
using std::vector;

DiscretizationPolyhedralFEMOptimized::
DiscretizationPolyhedralFEMOptimized(mesh::Mesh & grid,
                                     const FiniteElementConfig & config,
                                     const std::vector<int> & fracture_face_markers,
                                     const std::vector<int> & neumann_face_markers)
    : DiscretizationPolyhedralFEM::DiscretizationPolyhedralFEM(grid, config,
                                                               fracture_face_markers,
                                                               neumann_face_markers)

{}

void DiscretizationPolyhedralFEMOptimized::build_(mesh::Cell & cell)
{
  FiniteElementData const * p_master = nullptr;
  std::vector<size_t> order;
  if ( known_element_(cell, order, p_master) )
  {
    cell.reorder_vertices(order);
    scale_cell_fem_data_(cell, *p_master);
    throw std::runtime_error("Awesome");
  }
  else DiscretizationPolyhedralFEM::build_(cell);
}

bool DiscretizationPolyhedralFEMOptimized::known_element_(mesh::Cell const & cell,
                                                          std::vector<size_t> & order,
                                                          FiniteElementData const *& p_master) const
{
  size_t const hsh = cell.n_vertices();
  if (_cell_data_compressed.count(hsh))
  {
    auto it = _cell_data_compressed.find(hsh);
    for (auto const & master : it->second)
    {
      auto result = Isomorphism::check(*master.topology, *cell.polyhedron());
      if (result.first)
      {
        order = result.second;
        p_master = &master;
        return true;
      }
    }
  }
  return false;
}

void DiscretizationPolyhedralFEMOptimized::build_cell_data_(mesh::Cell const & cell)
{
  DiscretizationFEMBase::build_cell_data_(cell);
  // remember topology
  auto & it = _cell_data_compressed[cell.n_vertices()];
  it.emplace_back();
  auto & saved = it.back();
  saved = _cell_data[cell.index()];
  saved.topology = cell.polyhedron();
}

void DiscretizationPolyhedralFEMOptimized::scale_cell_fem_data_(mesh::Cell const & cell,
                                                                FiniteElementData const & master)
{
  size_t const nv = cell.n_vertices();
  size_t const nq = master.points.size();
  FiniteElementData data(nv, nq);

  angem::Tensor2<3, double> du_dx;
  for (size_t q = 0; q < nq; ++q)
  {
    double detJ;
    for (auto grad : master.points[q].grads)
      std::cout <<  grad << std::endl;
    compute_detJ_and_invert_cell_jacobian_(master.points[q].grads, du_dx, detJ,
                                           cell.vertex_coordinates());
    if (detJ <= 0) throw std::runtime_error("Transformation det(J) is negative " + std::to_string(detJ));
    update_shape_grads_(master.points[q].grads, du_dx, data.points[q].grads);
  }
}

void DiscretizationPolyhedralFEMOptimized::
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

void DiscretizationPolyhedralFEMOptimized::
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


}  // end namespace discretization
