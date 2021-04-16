#include "DiscretizationPolyhedralFEMOptimized.hpp"
#include "Elements/PolyhedralElementScaled.hpp"
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
  std::vector<size_t> order;  // vertex ordering in case we need to reorder
  if ( auto master = known_element_(cell, order) )
  {
    std::cout << "found similar element " <<  cell.index() << std::endl;
    cell.reorder_vertices(order);
    _element = std::make_shared<PolyhedralElementScaled>( cell, _grid, *master, _config );
  }
  else {
    std::cout << "new element topology " <<  cell.index() << std::endl;
    DiscretizationPolyhedralFEM::build_(cell);
    // remember topology for future re-use
    _masters[cell.n_vertices()].push_back(
        std::dynamic_pointer_cast<PolyhedralElementBase>(_element ));
 }
}

std::shared_ptr<PolyhedralElementBase>
DiscretizationPolyhedralFEMOptimized::known_element_(mesh::Cell const & cell,
                                                     std::vector<size_t> & order) const
{
  size_t const hsh = cell.n_vertices();
  if (_masters.count(hsh))
  {
    auto it = _masters.find(hsh);
    for (auto master : it->second)
    {
      auto result = Isomorphism::check(*master->host_topology(),
                                       *cell.polyhedron());
      if (result.first)
      {
        order = result.second;
        return master;
      }
    }
  }
  std::cout << "nada" << std::endl;
  return nullptr;
}

}  // end namespace discretization
