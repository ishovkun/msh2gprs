#include "PolyhedralElementScaled.hpp"
#include "logger/Logger.hpp"

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
  compute_detJ_and_invert_cell_jacobian_(master.grads, du_dx, detJ, vertex_coord);

  target.weight = detJ * master.weight;
  target.values = master.values;

  if (detJ <= 0)
    throw std::runtime_error("Cell Transformation det(J) is negative " +
                             std::to_string(detJ));
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
   build_fe_face_point_data_(vert_coord, master_data.center,
                             face_data.center, du_dx, basis);
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

  double detJ;
  compute_detJ_and_invert_face_jacobian_(master.grads, du_dx, detJ, vertex_coord, basis);
  if ( detJ <= 0 ) throw std::runtime_error("Face Transformation det(J) is negative " + std::to_string(detJ));
  update_shape_grads_(master.grads, du_dx, current.grads);
}

void PolyhedralElementScaled::
compute_detJ_and_invert_face_jacobian_(const std::vector<angem::Point<3,double>> & ref_grad,
                                       angem::Tensor2<3, double> & J_inv,
                                       double & detJ,
                                       std::vector<angem::Point<3,double>> const & vertex_coord,
                                       const angem::Basis<3,double> & basis) const
{
  angem::Plane<double> plane(vertex_coord);
  // {
  //   std::cout << "basis.norm = " << basis(2)  << std::endl;
  //   std::cout << "face.norm  = " << plane.normal()  << std::endl;
  // }

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
  detJ = det(J);
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
    if (it == master_verts.end())
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

}  // end namespace discretization
