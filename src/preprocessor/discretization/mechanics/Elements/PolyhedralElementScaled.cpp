#include "PolyhedralElementScaled.hpp"
#include "logger/Logger.hpp"
#include "PFEM_integration/IntegrationRule3d.hpp"
#include "PFEM_integration/IntegrationRule2d.hpp"
#include "PFEM_integration/IntegrationRuleFracture.hpp"

namespace discretization {

using Point = angem::Point<3,double>;

#ifdef WITH_EIGEN
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
  auto const vert_coord = _parent_cell.vertex_coordinates();
  auto const & rule = _master.integration_rule3();
  _cell_data = rule.integrate(_parent_cell.vertex_coordinates());
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

FiniteElementData PolyhedralElementScaled::get_face_data(size_t iface)
{
  auto const * const face = _parent_cell.faces()[iface];
  auto const * const master_face = _master.host_cell().faces()[iface];
  auto vert_coord = _parent_cell.faces()[iface]->vertex_coordinates();
  reorder_face_vertices_(*face, *master_face, vert_coord);

  auto basis = get_face_basis_(*face, host_cell());
  auto const & rule = _master.integration_rule2(iface);

  return rule.integrate( vert_coord, basis );

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
}
FiniteElementData PolyhedralElementScaled::get_fracture_data(const size_t iface,
                                                             const angem::Basis<3,double> basis)
{
  return _master.integration_rule_frac(iface).integrate(_parent_cell.vertex_coordinates(), basis);
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

#endif
}  // end namespace discretization
