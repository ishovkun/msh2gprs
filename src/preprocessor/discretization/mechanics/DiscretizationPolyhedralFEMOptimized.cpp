#include "DiscretizationPolyhedralFEMOptimized.hpp"
#include "Isomorphism.hpp"         // provides isomorphism


namespace discretization {

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
  FiniteElementDataTopology const * p_master = nullptr;
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
                                                          FiniteElementDataTopology const *& p_master) const
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
    return false;
  }
  else return false;
}

void DiscretizationPolyhedralFEMOptimized::scale_cell_fem_data_(mesh::Cell const & cell,
                                                                FiniteElementDataTopology const & data)
{

}


}  // end namespace discretization
