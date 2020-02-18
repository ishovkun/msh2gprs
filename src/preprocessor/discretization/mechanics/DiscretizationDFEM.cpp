#include "DiscretizationDFEM.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include <stdexcept>
#ifdef WITH_EIGEN
#include "DFEMElement.hpp"
#endif

namespace discretization
{
using Point = angem::Point<3,double>;
using api = gprs_data::GmshInterface;

DiscretizationDFEM::DiscretizationDFEM(const mesh::Mesh &grid)
    : _grid(grid)
{
  #ifndef WITH_GMSH
  throw std::runtime_error("Cannot use DFEM method without linking to GMsh");
  #endif
  #ifndef WITH_EIGEN
  throw std::runtime_error("Cannot use DFEM method without linking to Eigen");
  #endif
}

#ifdef WITH_GMSH

void DiscretizationDFEM::build()
{
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
  {
    api::initialize_gmsh();
    std::cout << "cell.index() = " << cell->index() << std::endl;
    DFEMElement discr_element(*cell);
    api::finalize_gmsh();
    exit(0);
  }
}

#else
void DiscretizationDFEM::build() {}
#endif

}  // end namepsace discretization
