#include "DiscretizationDFEM.hpp"

namespace gprs_data
{
using Point = angem::Point<3,double>;

DiscretizationDFEM::DiscretizationDFEM(const mesh::Mesh & grid)
    : _grid(grid)
{}

#ifdef WITH_GMSH

void DiscretizationDFEM::build()
{
  // build_gmsh_simple();
  // exit(0);
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
  {
    std::cout << "cell.index() = " << cell->index() << std::endl;
    build_(*cell);
    exit(0);
  }
}

void DiscretizationDFEM::build_(const mesh::Cell & cell)
{
  GmshInterface::initialize_gmsh();
  GmshInterface::build_triangulation(cell);
  build_shape_functions_();
  GmshInterface::finalize_gmsh();
}

void DiscretizationDFEM::build_shape_functions_()
{
  std::vector<int> element_types;
  std::vector<std::vector<std::size_t> > element_tags;
  std::vector<std::vector<std::size_t> > node_tags;
  gmsh::model::mesh::getElements(element_types, element_tags, node_tags);
  // std::vector<int> volume_types = get_msh_volume_types_(element_types);
  std::cout << "tags" << std::endl;
  for (auto type : element_types)
    std::cout << type << std::endl;

  // 1. make sure that we only have tetras (look at element_types)
  // 2. ask gmsh to provide gaussian points, shape functions, and jacobians
  // for our tetras
  // 3. build element jacobian for the homogeneous laplace equation
}

#else
void DiscretizationDFEM::build() {}
#endif

}  // end namepsace gprs_data
