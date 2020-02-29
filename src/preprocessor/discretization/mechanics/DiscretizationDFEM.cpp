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
  int cnt = 0;
  {
    auto & cell = _grid.cell(100);
    // for (auto v : cell.vertices())
    //   std::cout << _grid.vertex(v) << "\t|\t";
    // std::cout << std::endl;
    // auto c = cell.polyhedron();
    // for (auto v : c->get_points())
    //   std::cout << v << "\t|\t";
    // std::cout << std::endl;
    // exit(0);
    // std::cout << "c->id() = " << c->id() << std::endl;
    api::initialize_gmsh();
    DFEMElement discr_element(cell, _msrsb_tol);
    mesh::Mesh _element_grid;
    api::finalize_gmsh();
    exit(0);
    // if (cnt++ > 2) break;
  }

  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
  {
    std::cout << "cell->index() = " << cell->index() << std::endl;
    api::initialize_gmsh();
    DFEMElement discr_element(*cell, _msrsb_tol);
    mesh::Mesh _element_grid;
    api::finalize_gmsh();
    // // exit(0);
    // if (cnt++ > 2) break;
  }
}

#else
void DiscretizationDFEM::build() {}
#endif

}  // end namepsace discretization
