#include "DiscretizationDFEM.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include <stdexcept>
#ifdef WITH_EIGEN
#include "DFEMElement.hpp"
#include "PolyhedralElementDirect.hpp"
#endif

namespace discretization
{
using Point = angem::Point<3,double>;
using api = gprs_data::GmshInterface;

DiscretizationDFEM::DiscretizationDFEM(const mesh::Mesh & grid, const double msrsb_tol)
    : _grid(grid), _msrsb_tol( msrsb_tol )
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
  // {
  //   auto & cell = _grid.cell(100);

  //   api::initialize_gmsh();
  //   PolyhedralElementDirect discr_elemnelement(cell);

  //   // DFEMElement discr_element(cell, _msrsb_tol);
  //   // mesh::Mesh _element_grid;
  //   api::finalize_gmsh();
  //   exit(0);
  //   // if (cnt++ > 2) break;
  // }

  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
  {
    std::cout << "cell->index() = " << cell->index() << std::endl;
    api::initialize_gmsh();
    /* DFEMElement discr_element(*cell, _msrsb_tol); */
    PolyhedralElementDirect discr_element(*cell);
    /* mesh::Mesh _element_grid; */
    api::finalize_gmsh();

   FiniteElementData cell_fem_data =  discr_element.get_cell_data();
   cell_fem_data.element_index = cell->index();
   _cell_data.push_back( cell_fem_data );

    // // exit(0);
    // if (cnt++ > 2) break;
  }
}

#else
void DiscretizationDFEM::build() {}
#endif

}  // end namepsace discretization
