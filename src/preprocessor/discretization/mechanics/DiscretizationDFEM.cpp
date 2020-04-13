#include "DiscretizationDFEM.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include <stdexcept>
#ifdef WITH_EIGEN
#include "DFEMElement.hpp"
#include "PolyhedralElementDirect.hpp"
#endif
#include "StandardFiniteElement.hpp"
#include "ProgressBar.hpp"  // provides ProgressBar


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
  {
    auto & cell = _grid.cell(0);

    api::initialize_gmsh();
    PolyhedralElementDirect de(cell);
    StandardFiniteElement fe(cell);
    // de.debug_save_shape_functions_("cell20.vtk");
    // de.debug_save_boundary_face_solution("cell20_faces.vtk");
    // gmsh::write("cell20.msh");

    // DFEMElement discr_element(cell, _msrsb_tol);
    // mesh::Mesh _element_grid;
    api::finalize_gmsh();
    exit(0);
  }

  _face_data.resize( _grid.n_faces() );

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

   std::vector<FiniteElementData> face_data = discr_element.get_face_data();
   size_t iface = 0;
   for ( const mesh::Face * face : cell->faces() )
   {
     if ( _face_data[face->index()].points.empty() )
     {
       face_data[iface].element_index = face->index();
       _face_data[face->index()] = face_data[iface++];
     }
   }

    // // exit(0);
    // if (cnt++ > 2) break;
  }
}

#else
void DiscretizationDFEM::build() {}
#endif

}  // end namepsace discretization
