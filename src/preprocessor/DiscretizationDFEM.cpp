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
  build_jacobian_();
  // run_msrsb_process();
  // extract_shape_functions_();

}

void DiscretizationDFEM::build_jacobian_()
{
  // 1. ask gmsh to provide gaussian points, shape functions, and jacobians
  // for our tetras
  // 2. build element jacobian for the homogeneous laplace equation
  std::vector<int> element_types;
  std::vector<std::vector<std::size_t> > element_tags;
  std::vector<std::vector<std::size_t> > node_tags;
  GmshInterface::get_elements(element_types, element_tags, node_tags, /* dim = */ 3);
  std::cout << "element_tags.size() = " << element_tags.size() << std::endl;
  std::cout << "node-tags.size() = " << node_tags.size() << std::endl;
  for (std::size_t itype = 0; itype < element_types.size(); ++itype)
  {
    const size_t type = element_types[itype];
    for (const size_t tag : element_tags[itype])
    {
      build_local_matrix_(tag);
    }
  }
}

void DiscretizationDFEM::build_local_matrix_(const size_t element_tag)
{
  
}

#else
void DiscretizationDFEM::build() {}
#endif

}  // end namepsace gprs_data
