#include "DiscretizationDFEM.hpp"
#include "gmsh_interface/FeValues.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

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

  Eigen::SparseMatrix<double,Eigen::RowMajor> mat(node_tags.size(), node_tags.size());
  for (std::size_t itype = 0; itype < element_types.size(); ++itype)
  {
    const size_t type = element_types[itype];
    FeValues fe_values(type, element_tags[itype].size());

    for (const size_t tag : element_tags[itype])
    {
      fe_values.update(tag);
      Eigen::MatrixXd cell_matrix = Eigen::MatrixXd::Zero(fe_values.n_vertices(),
                                                          fe_values.n_vertices());
      std::cout << "cell_matrix = " << cell_matrix << std::endl;

      for (size_t q = 0; q < fe_values.n_q_points(); ++q)
      {
        for (size_t i = 0; i < fe_values.n_vertices(); ++i)
          for (size_t j = 0; j < fe_values.n_vertices(); ++j)
            cell_matrix(i, j) += (fe_values.grad(i, q) * // grad phi_i(x_q)
                                  fe_values.grad(j, q) * // grad phi_j(x_q)
                                  fe_values.JxW(q));           // dx
      }
      // distribute_local_to_global_();
      std::cout << "cell_matrix = " << cell_matrix << std::endl;
      exit(0);
    //   build_local_matrix_(type, tag);
    }
  }
}

#include "gmsh.h"

void DiscretizationDFEM::build_local_matrix_(const int element_type,
                                             const size_t element_tag)
{
  std::cout << "element_tag = " << element_tag << std::endl;
  std::cout << "element_type = " << element_type << std::endl;
  exit(0);
}

#else
void DiscretizationDFEM::build() {}
#endif

}  // end namepsace gprs_data
