#pragma once
#include "angem/VTK_ID.hpp"
#include "mesh/Mesh.hpp"
#include "angem/Tensor2.hpp"

namespace discretization {

using angem::VTK_ID;
using Point = angem::Point<3,double>;
template<VTK_ID id> constexpr size_t N_ELEMENT_VERTICES;

/**
 * This is an abstract class that combines methods for standard finite elements.
 * It can be used to compute the values of shape functions, their gradients,
 * and JxW values.
 */
template<VTK_ID vtk_id>
class FeValues
{
 public:
  /**
   * Constructor.
   * This is an abstract class so the constuctor is protected.
   */
  FeValues(const mesh::Mesh & grid);
  /**
   * Update the needed quantities in the cell.
   * Use this method before using the fem quantities in the new cell during the loop.
   */
  void update(const mesh::Cell & cell);
  /**
   * Update the needed quantities in the cell in the gived integration point.
   * The point coordinate is the real (not reference) coordinates of the integration point.
   */
  void update(const mesh::Cell & cell, const angem::Point<3,double> & point);
  /**
   * Update the needed quantities in the face.
   * Use this method before using the fem quantities in the new face during the loop.
   */
  void update(const mesh::Face & face);
  /**
   * Update the needed quantities in the face in the gived integration point.
   * The point coordinate is the real (not reference) coordinates of the integration point.
   */
  void update(const mesh::Face & cell, const angem::Point<3,double> & point);
  /**
   * Computes a vector of integration points in the master element
   */
  std::vector<Point> get_master_integration_points() const;
  /**
   * Computes a vector of integration point weights in the master element
   */
  std::vector<double> get_master_integration_weights() const;

  /**
   * Return the value of the shape function indexed by shape_index in the
   * qpoint integration point in the current element.
   * NOTE: must call update() before calling this function.
   */
  double value(const size_t shape_index, const size_t qpoint) const;
  /**
   * Return the gradient of the shape function indexed by shape_index in the
   * qpoint integration point in the current element.
   * NOTE: must call update() before calling this function.
   */
  Point grad(const size_t shape_index, const size_t qpoint) const;
  /**
   * Return the JxW value at the integration point indexed by qpoint.
   * JxW is a product of the determinant of the transformation jacobian by
   * the weight of the integration point.
   * NOTE: must call update() before calling this function.
   */
  double JxW(const size_t qpoint) const;

 protected:
  /**
   * Given element vertices and integration points, update shape values,
   * gradients, and the determinant of the trasformation jacobian
   */
   void update_();
  /**
   * Evalues the values of the ith shape function in point in the reference element
   * \param[in] point  : coordinates in the reference element
   * \param[in] vertex : the local number of the shape function
   * returns the value of the given shape function in a given point
   */
  double eval_(const Point & point, const size_t vertex) const;
  /**
   * Evalues the derivatives of the ith shape function in point in the reference element
   * \param[in] point  : coordinates in the reference element
   * \param[in] vertex : the local number of the shape function
   * returns the derivative of the given shape function in a given point
   */
   Point eval_derivative_(const Point & point, const size_t vertex) const;

  /**
   * Compute transformation jacobian matrix dx / dx_ref
   * Input:
   * \param[in] ref_grad : shape function grad in current q-point in master element
   * Returns transposed matrix dx / dx_ref
   */
   angem::Tensor2<3, double> compute_jacobian_(const std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & ref_grad) const;

  // compute shape values in q-points in master element
  void update_shape_values_();
  // compute the gradients of shape fucntions in qpoint in master element
  std::array<Point,N_ELEMENT_VERTICES<vtk_id>> compute_ref_gradient_(const size_t qpoint) const;

  /**
   * Map from real coordinates xyz to local coordinates in the master element uvw
   */
   Point map_real_to_local_(const Point & p) const;

  /**
   * Get vertex coordinate
   */
  const Point & vertex_(const size_t i) const {return _grid.vertex(_element_vertices[i]);}

  const mesh::Mesh & _grid;
  std::vector<Point> _qpoints;
  std::vector<std::array<double,N_ELEMENT_VERTICES<vtk_id>>> _shape_values;
  std::vector<std::array<Point,N_ELEMENT_VERTICES<vtk_id>>> _shape_grads;
  std::vector<double> _weights;
  std::vector<double> _determinants;
  std::vector<size_t> _element_vertices;
};

template<VTK_ID vtk_id>
FeValues<vtk_id>::FeValues(const mesh::Mesh & grid)
    : _grid(grid)
{}


template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(const mesh::Cell & cell)
{
  static_assert(vtk_id == VTK_ID::TetrahedronID,
                "This function only exists for 3D elements and is only implemented for tetras");
  _element_vertices = cell.vertices();
  _qpoints = get_master_integration_points();
  _weights = get_master_integration_weights();
  update_();
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(const mesh::Cell & cell, const angem::Point<3,double> & point)
{
  static_assert(vtk_id == VTK_ID::TetrahedronID,
                "This function only exists for 3D elements and is only implemented for tetras");
  _element_vertices = cell.vertices();
  _weights = {1.0};
  _qpoints.resize(1);
  _qpoints[0] = point;
  assert(false && "Implement mapping of real coordinates into reference coordinated");
  update_();
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(const mesh::Face & face)
{
  static_assert(vtk_id == VTK_ID::TriangleID,
                "This function only exists for 2D elements and is only implemented for triangles");

  _element_vertices = face.vertices();
  _qpoints = get_master_integration_points();
  _weights = get_master_integration_weights();
  update_();
}

template<VTK_ID vtk_id>
angem::Tensor2<3, double> FeValues<vtk_id>::compute_jacobian_(const std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & ref_grad) const
{
  // compute jacobian
  // see Becker E., Carey G., Oden J. Finite elements. An Introduction Volume 1 1981
  // Eq. 5.3.6
  // NOTE: this is a transposed jacobian matrix
  angem::Tensor2<3, double> J;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      for (std::size_t v=0; v<N_ELEMENT_VERTICES<vtk_id>; ++v)
        J( j, i ) += ref_grad[v][j] * _grid.vertex(_element_vertices[v])[i];
  return J;
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update_()
{
  // compute shape values in q-points in master element
  update_shape_values_();

  _determinants.assign(_qpoints.size(), 0.0);
  _shape_grads.assign(_qpoints.size(), std::array<Point,N_ELEMENT_VERTICES<vtk_id>>());
  // loop integration points
  for (std::size_t q=0; q<_qpoints.size(); ++q)
  {
    // compute shape grad in in master element
    const auto ref_grad = compute_ref_gradient_(q);
    // compute transformation jacobian
    angem::Tensor2<3, double> dx_du = compute_jacobian_(ref_grad);
    // compute the determinant of transformation jacobian
    _determinants[q] = det(dx_du);
    // must be positive
    if ( _determinants[q] <= 0 )
      throw std::runtime_error("Transformation Jacobian is not invertible");
    // invert the jacobian to compute shape function gradients
    angem::Tensor2<3, double> du_dx = invert(dx_du);

    // compute the true shape function gradients
    for (size_t vertex = 0; vertex < N_ELEMENT_VERTICES<vtk_id>; ++vertex)
      for (size_t i=0; i<3; ++i)
        for (size_t j = 0; j < 3; ++j)
          // d phi_vert / dx_i = (d phi_vert / d u_j) * (d u_j / d x_i)
          _shape_grads[q][vertex][i] += ref_grad[vertex][j] * du_dx(j, i);
  }

}

template<VTK_ID vtk_id>
double FeValues<vtk_id>::value(const size_t shape_index, const size_t qpoint) const
{
  assert( qpoint < _qpoints.size() && "qpoint too large" );
  assert( shape_index < N_ELEMENT_VERTICES<vtk_id> && "shape_index too large");
  return _shape_values[qpoint][shape_index];
}

template<VTK_ID vtk_id>
Point FeValues<vtk_id>::grad(const size_t shape_index, const size_t qpoint) const
{
  assert( qpoint < _qpoints.size() && "qpoint too large" );
  assert( shape_index < N_ELEMENT_VERTICES<vtk_id> && "shape_index too large");
  return _shape_grads[qpoint][shape_index];
}

template<VTK_ID vtk_id>
double FeValues<vtk_id>::JxW(const size_t qpoint) const
{
  assert( qpoint < _qpoints.size() && "qpoint too large" );
  return _determinants[qpoint] * _weights[qpoint];
}


template<VTK_ID vtk_id>
void FeValues<vtk_id>::update_shape_values_()
{
  _shape_values.resize( _qpoints.size() );
  for (std::size_t q=0; q<_qpoints.size(); ++q)
    for (std::size_t v=0; v<N_ELEMENT_VERTICES<vtk_id>; ++v)
      _shape_values[q][v] = eval_(_qpoints[q], v);
}

template<VTK_ID vtk_id>
std::array<Point,N_ELEMENT_VERTICES<vtk_id>> FeValues<vtk_id>::compute_ref_gradient_(const size_t q) const
{
  std::array<Point,N_ELEMENT_VERTICES<vtk_id>> ref_grad;
  for (std::size_t v=0; v<N_ELEMENT_VERTICES<vtk_id>; ++v)
    ref_grad[v] = eval_derivative_(_qpoints[q], v);
  return ref_grad;
}


}  // end namespace discretization

#include "FeValuesTriangle.hpp"
#include "FeValuesTetra.hpp"
