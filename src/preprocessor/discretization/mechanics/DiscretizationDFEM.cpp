#include "DiscretizationDFEM.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#ifdef WITH_EIGEN
#include "DFEMElement.hpp"
#endif

namespace discretization
{
using Point = angem::Point<3,double>;
using api = gprs_data::GmshInterface;

DiscretizationDFEM::DiscretizationDFEM(const mesh::Mesh & grid)
    : _grid(grid)
{}

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
