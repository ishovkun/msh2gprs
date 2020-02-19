#include "FeValues.hpp"
#include "../utils.hpp"

#ifdef WITH_GMSH
#include "gmsh.h"
#endif

#include <cassert>
#include <stdexcept>
#include <iostream>

namespace gprs_data {

FeValues::FeValues(const int element_type, const size_t n_elements)
    : _n_comp(1), _element_type(element_type), _n_elements(n_elements)
{
  initialize_();
}

#ifdef WITH_GMSH
void FeValues::initialize_()
{
  gmsh::model::mesh::getIntegrationPoints(_element_type, "Gauss1", _ref_points, _weights);
  gmsh::model::mesh::getBasisFunctions(_element_type, _ref_points,
                                       "GradLagrange", _n_comp, _ref_gradients);
  gmsh::model::mesh::getJacobians(_element_type, _ref_points,
                                  _jacobians, _determinants, _true_points,
                                  /* tag = */ -1, /* task = */ 0, /* n_tasks = */ 1);
}

#else
void FeValues::initialize_()
{}
#endif

void FeValues::update(const size_t cell_number)
{
  const size_t dim = 3;
  const double nv = n_vertices();
  _grad.assign( _ref_gradients.size() * n_q_points(), 0.0 );
  _inv_determinants.resize(_weights.size());
 
  for (std::size_t q=0; q<n_q_points(); ++q)
  {
    // build local dx_du matrix
    // dim*dim derivative entries for n_gauss_points for each cell
    const size_t j_beg = ( cell_number * n_q_points() + q  ) * (dim * dim);
    const size_t j_end = ( cell_number * n_q_points() + q+1) * (dim * dim);
    std::vector<double> dx_du (_jacobians.begin() + j_beg, _jacobians.begin() + j_end);
    // gmsh stores jacobians in transposed format
    // dx_du = transpose3x3(dx_du);

    // invert du_dx = inv(dx_du)
    std::vector du_dx(dim*dim, 0.0);
    invert_matrix_3x3(dx_du, du_dx);
    _inv_determinants[q] = determinant_3x3(du_dx);

    // compute shape function gradients
    for (size_t vertex = 0; vertex < nv; ++vertex)
      for (std::size_t i=0; i<dim; ++i)
        for (std::size_t j=0; j<dim; ++j)
          _grad[q*dim*nv + dim*vertex + i] +=
              du_dx[i*dim + j] * _ref_gradients[q*dim*nv + dim*vertex + j];
  }
  // debug_print_cell_config();
}

double FeValues::value(size_t i, size_t q) const
{
  throw std::runtime_error("not implemeneted");
  const double dim = 3;
  const double n_vertices = _ref_gradients.size() / dim / n_q_points();
  assert( q < n_q_points() );

  return 0.0;
}

size_t FeValues::n_q_points() const
{
  return _weights.size();
}

size_t FeValues::n_vertices() const
{
  const double dim = 3;
  return _ref_gradients.size() / dim / n_q_points();
}


Point FeValues::grad(size_t vertex, size_t q) const
{
  const double dim = 3;
  const double nv = n_vertices();
  assert( q < n_q_points() );
  assert( vertex < nv );
  Point p;
  p[0] = _grad[q*dim*nv + dim*vertex];
  p[1] = _grad[q*dim*nv + dim*vertex + 1];
  p[2] = _grad[q*dim*nv + dim*vertex + 2];
  return p;
}

double FeValues::JxW(const size_t q) const
{
  // return _weights[q] * _inv_determinants[q];
  return _weights[q] * _determinants[q];
}

void FeValues::debug_print_cell_config()
{
  std::cout << "printing cell" << std::endl;
  size_t cell_number = 0;
  const size_t j_beg = cell_number       * 1 * (9);
  const size_t j_end = (cell_number + 1) * 1 * (9);
  std::vector<double> dx_du (_jacobians.begin() + j_beg, _jacobians.begin() + j_end);
  std::cout << "direct jac" << std::endl;
  for (auto v : dx_du)
    std::cout << v << " ";
  std::cout << std::endl;
  dx_du = transpose3x3(dx_du);
  std::cout << "direct transpose jac" << std::endl;
  for (auto v : dx_du)
    std::cout << v << " ";
  std::cout << std::endl;
  std::cout << "integration points" << std::endl;
  for (auto p : _ref_points)
    std::cout << p << " ";
  std::cout << std::endl;

  // gmsh::model::mesh::getJacobians(_element_type, _ref_points,
  //                                 _jacobians, _determinants, _true_points,
  //                                 /* tag = */ -1, /* task = */ 0, /* n_tasks = */ 1);


  exit(0);
}

}  // end namespace gprs_data
