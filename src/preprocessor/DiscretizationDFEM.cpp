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
      build_local_matrix_(type, tag);
    }
  }
}

#include "gmsh.h"

void DiscretizationDFEM::build_local_matrix_(const int element_type,
                                             const size_t element_tag)
{
  std::cout << "element_tag = " << element_tag << std::endl;
  std::cout << "element_type = " << element_type << std::endl;
  std::vector<double> integration_points;
  std::vector<double> integration_weights;
  gmsh::model::mesh::getIntegrationPoints(element_type, "Gauss1", integration_points, integration_weights);

  // std::cout << "integration_points.size() = " << integration_points.size() << std::endl;
  // for (std::size_t i=0; i<integration_points.size()/3; ++i)
  // {
  //   std::cout << integration_points[3*i] << " "
  //             << integration_points[3*i+1] << " "
  //             << integration_points[3*i+2] << " "
  //             << integration_weights[i]
  //             << std::endl;
  // }

  // const double n_vertices = GmshInterface::get_n_vertices(element_type);
  std::vector<double> grad_phi;  // basis function gradients

  // this will be GmshInterface::getFunctionGradients()
  int n_comp = 1;
  gmsh::model::mesh::getBasisFunctions(element_type, integration_points, "GradLagrange",
                                       n_comp, grad_phi) ;

  std::cout <<  " grad u" << std::endl;
  for (auto v : grad_phi)
    std::cout << v << std::endl;

  // c) getBasisFunctions (or getBasisFunctionsForElements)
  // d) getJacobians probably contains both sf gradients and JxW values
  std::vector<double> jacobians, determinants, points;
  gmsh::model::mesh::getJacobians(element_type,
                                 integration_points,
                                 jacobians,
                                 determinants,
                                 points,
                                 /* tag = */ -1,
                                 /* task = */ 0,
                                 /* n_tasks = */ 1);
  //
        // vShapeDerivatives[irow * points_in_element + inodes] +=
        //   vInversedJac[irow * stdelement[element_form].dimension_size + icol] * (*vPointDeriv_)[icol];
  exit(0);
}

#else
void DiscretizationDFEM::build() {}
#endif

}  // end namepsace gprs_data
