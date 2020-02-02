#include "DiscretizationDFEM.hpp"
#include "DFEMElement.hpp"

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
    exit(0);
    api::finalize_gmsh();
  }
}

#else
void DiscretizationDFEM::build() {}
#endif

}  // end namepsace discretization
