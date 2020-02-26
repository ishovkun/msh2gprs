#include "FeValues.hpp"
#include "GmshInterface.hpp"
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
  _cell_index = std::numeric_limits<size_t>::max();
  /* NOTE: this updates quantities in all the grid elements */
  gmsh::model::mesh::getIntegrationPoints(_element_type, "Gauss1", _ref_points, _weights);
  gmsh::model::mesh::getBasisFunctions(_element_type, _ref_points, "Lagrange", _n_comp, _ref_values);
  gmsh::model::mesh::getBasisFunctions(_element_type, _ref_points, "GradLagrange", _n_comp, _ref_gradients);
  gmsh::model::mesh::getJacobians(_element_type, _ref_points, _jacobians, _determinants, _true_points,
                                  /* tag = */ -1, /* task = */ 0, /* n_tasks = */ 1);
  get_elements_();
}

void FeValues::get_elements_()
{
  std::vector<int> element_types;
  std::vector<std::vector<size_t>> element_node_tags;
  std::vector<std::vector<size_t>> tags;
  GmshInterface::get_elements(element_types, tags, element_node_tags, 3, -1);
  _element_tags = tags[0];
}

#else
void FeValues::initialize_()
{}
#endif

void FeValues::update(const size_t cell_number)
{
  const size_t dim = 3;
  const double nv = n_vertices();
  _cell_index = cell_number;
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
}

void FeValues::update(const size_t cell_number, const std::vector<Point> & points)
{
  _cell_index = cell_number;
  const size_t tag = _element_tags[_cell_index];
  std::vector<double> local_points;
  local_points.reserve( 3*points.size() );
  for (const auto & p : points)
  {
    double u, v, w;
    gmsh::model::mesh::getLocalCoordinatesInElement(tag, p(0), p(1), p(2), u, v, w);
    local_points.push_back( u );
    local_points.push_back( v );
    local_points.push_back( w );
  }
  gmsh::model::mesh::getBasisFunctions(_element_type, local_points, "Lagrange",
                                       _n_comp, _ref_values);
  gmsh::model::mesh::getBasisFunctions(_element_type, local_points, "GradLagrange",
                                       _n_comp, _ref_gradients);
  _true_points.clear();
  std::vector<double> loc_jac, loc_det;
  gmsh::model::mesh::getJacobians(_element_type, local_points, _jacobians, _determinants, _true_points,
                                  /* tag = */ -1, /* task = 0 */ _cell_index, /* n_tasks = */ _n_elements);
  update(cell_number);
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


Point FeValues::grad(const size_t vertex, const size_t q) const
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

double FeValues::value(const size_t vertex, const size_t q) const
{
  const double nv = n_vertices();
  assert( vertex < nv );
  assert( q < n_q_points() );
  return _ref_values[ q*nv + vertex ];
}

double FeValues::JxW(const size_t q) const
{
  // return _weights[q] * _inv_determinants[q];
  return _weights[q] * _determinants[q];
}

}  // end namespace gprs_data
