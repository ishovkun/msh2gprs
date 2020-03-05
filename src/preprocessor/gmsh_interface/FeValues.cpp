#include "FeValues.hpp"
#include "GmshInterface.hpp"
#include "angem/Tensor2.hpp"

#include "../utils.hpp"

#ifdef WITH_GMSH
#include "gmsh.h"
#endif

#include <cassert>
#include <stdexcept>
#include <iostream>

namespace gprs_data {

FeValues::FeValues(const int element_type, const int element_tag)
    : _n_comp(1), _element_type(element_type), _element_tag(element_tag)
{
  initialize_();
}

#ifdef WITH_GMSH
void FeValues::initialize_()
{
  _cell_index = std::numeric_limits<size_t>::max();
  /* NOTE: this updates quantities in all the grid elements */
  gmsh::model::mesh::getIntegrationPoints(_element_type, "Gauss1", _ref_points, _weights);
  // std::cout << "n_weights = " << _weights.size() << std::endl;
  // std::cout << "ref point " << _ref_points[0] << " " << _ref_points[1] << " " << _ref_points[2] << std::endl;
  // for (size_t i = 0; i < _ref_points.size() / 3; ++i) {
  //   _ref_points[i+0] = 0.58541;
  //   _ref_points[i+1] = 0.138197;
  //   _ref_points[i+2] = 0.138197;
  // }

  // std::cout << "ref point mod " << _ref_points[0] << " " << _ref_points[1]
  //           << " " << _ref_points[2] << std::endl;

  gmsh::model::mesh::getBasisFunctions(_element_type, _ref_points, "Lagrange", _n_comp, _ref_values);
  gmsh::model::mesh::getBasisFunctions(_element_type, _ref_points, "GradLagrange", _n_comp, _ref_gradients);
  gmsh::model::mesh::getJacobians(_element_type, _ref_points, _all_jacobians, _all_determinants, _true_points,
                                  /* tag = */ -1, /* task = */ 0, /* n_tasks = */ 1);

  get_elements_();
}

void FeValues::get_elements_()
{
  std::vector<size_t> element_node_tags;
  std::vector<size_t> tags;
  gmsh::model::mesh::getElementsByType( _element_type, _element_tags, element_node_tags, _element_tag );
  _n_elements = tags.size();
}

#else
void FeValues::initialize_()
{}
#endif

void FeValues::update(const size_t cell_number)
{
  const size_t dim = 3;
  const double nv = n_vertices();
  const size_t nq = n_q_points();
  _cell_index = cell_number;
  _grad.assign( _ref_gradients.size() * nq, 0.0 );
  _determinants.resize( nq, 0.0 );

  for (std::size_t q=0; q<nq; ++q)
  {
    // build local dx_du matrix
    // dim*dim derivative entries for n_gauss_points for each cell
    const size_t j_beg = 9 * ( cell_number * nq + q  );
    const size_t j_end = 9 * ( cell_number * nq + q + 1);
    _determinants[q] = _all_determinants[cell_number * nq + q];

    // take a portion of the jacobian vector that corresponds to q
    std::vector<double> jac_local(_all_jacobians.begin() + j_beg, _all_jacobians.begin() + j_end);
    angem::Tensor2<3, double> dx_du(jac_local);
    // invert jacobian
    angem::Tensor2<3, double> du_dx = invert(dx_du);

    for (size_t vertex = 0; vertex < nv; ++vertex)
      for (size_t i=0; i<dim; ++i)
      {
        const size_t offset = q * dim * nv;
        const size_t n = dim * vertex;
        for (size_t j = 0; j < dim; ++j)
          _grad[offset + n + i] += du_dx(i, j) * _ref_gradients[offset + n + j];
      }
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
  gmsh::model::mesh::getJacobians(_element_type, local_points, _all_jacobians, _determinants, _true_points,
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
  Point p (_grad[q*dim*nv + dim*vertex],
           _grad[q*dim*nv + dim*vertex + 1],
           _grad[q*dim*nv + dim*vertex + 2]);

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
  return _weights[q] * std::fabs(1/_determinants[q]);
  // return _weights[q] * std::fabs(_determinants[q]);
}

}  // end namespace gprs_data
