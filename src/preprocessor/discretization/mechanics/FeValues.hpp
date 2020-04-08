#pragma once

#include "mesh/Mesh.hpp"
#include "angem/VTK_ID.hpp"
#include "angem/Tensor2.hpp"

namespace discretization {

using angem::VTK_ID;
using Point = angem::Point<3,double>;
template<VTK_ID id> constexpr size_t N_ELEMENT_VERTICES;
template<VTK_ID id> constexpr size_t ELEMENT_DIM;

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
  FeValues();
  /**
   * Update the needed internal quantities in the cell.
   * Use this method before using the fem quantities in the new cell during the loop.
   */
  void update(const mesh::Cell & cell);
  /**
   * Update the needed internal quantities in the cell in the gived integration point.
   * The point coordinate is the real (not reference) coordinates of the integration point.
   */
  void update(const mesh::Cell & cell, const angem::Point<3,double> & point);
  /**
   * Update the needed internal quantities in the face element.
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
   * Return the value of the shape function indexed by shape_index in element center.
   * NOTE: must call update() before calling this function.
   */
  double value_center(const size_t shape_index) const;
  /**
   * Return the gradient of the shape function indexed by shape_index in the
   * qpoint integration point in the current element.
   * NOTE: must call update() before calling this function.
   */
  Point grad(const size_t shape_index, const size_t qpoint) const;
  /**
   * Return the gradient of the shape function indexed by shape_index in the
   * element center.
   * NOTE: must call update() before calling this function.
   */
  Point grad_center(const size_t shape_index) const;
  /**
   * Return the JxW value at the integration point indexed by qpoint.
   * JxW is a product of the determinant of the transformation jacobian by
   * the weight of the integration point.
   * NOTE: must call update() before calling this function.
   */
  double JxW(const size_t qpoint) const;
  /**
   * Return the detJ value at the element center
   * NOTE: must call update() before calling this function.
   */
  double detJ_center() const;

  /**
   * Number of integration points in the quadrature rule.
   */
  size_t n_integration_points() const { return _qpoints.size(); }

 protected:
  /**
   * Given element vertices and integration points, update shape values,
   * gradients, and the determinant of the trasformation jacobian.
   * This function is for cells (3D elements).
   */
  void update_();

  /**
   * Update data in a single point
   */
  void update_(const Point & p,
               std::array<double,N_ELEMENT_VERTICES<vtk_id>> & shape_values,
               std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & shape_grads,
               double & determinant) const;
  /**
   * Given a vector of the grid vertex indices, fill out the interal
   * _vertex_coord array.
   */
  void update_vertex_coord_(const std::vector<Point> & vertex_coordinates);
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
  void update_shape_values_(const Point                                   & p,
                            std::array<double,N_ELEMENT_VERTICES<vtk_id>> & values) const;
  // compute shape grads in q-points in master element
  void update_shape_grads_(const std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & ref_grad,
                           const angem::Tensor2<3, double> & du_dx,
                           std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & grads) const;
  // compute the gradients of shape fucntions in qpoint in master element
  std::array<Point,N_ELEMENT_VERTICES<vtk_id>> compute_ref_gradient_(const Point &p) const;

  /**
   * Map from real coordinates xyz to local coordinates in the master element uvw
   */
  Point map_real_to_local_(const Point & p) const;

  /**
   * Compute cell jacobian, its determinant, and invert it.
   */
  void compute_detJ_and_invert_cell_jacobian_(const std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & ref_grad,
                                              angem::Tensor2<3, double> & du_dx,
                                              double & detJ) const;
  /**
   * Compute face jacobian, its determinant, and invert it.
   * Need a separate method for it since the jacobian will be 2x2 instead of 3x3.
   */
  void compute_detJ_and_invert_face_jacobian_(const std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & ref_grad,
                                              angem::Tensor2<3, double> & du_dx,
                                              double & detJ) const;

  /*  ---------------------------------Variables ---------------------------------- */

  std::vector<Point> _qpoints;  // gauss point coordinates in master element
  Point              _center;   // center coordinates in master element
  std::vector<std::array<double,N_ELEMENT_VERTICES<vtk_id>>> _shape_values;
  std::vector<std::array<Point,N_ELEMENT_VERTICES<vtk_id>>> _shape_grads;
  std::array<double,N_ELEMENT_VERTICES<vtk_id>> _shape_values_center;
  std::array<Point,N_ELEMENT_VERTICES<vtk_id>>  _shape_grads_center;
  std::vector<double> _weights;
  std::vector<double> _determinants;
  double              _determinant_center;
  std::array<Point,N_ELEMENT_VERTICES<vtk_id>> _vertex_coord;
};

template<VTK_ID vtk_id>
FeValues<vtk_id>::FeValues()
{}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(const mesh::Cell & cell)
{
  static_assert(ELEMENT_DIM<vtk_id> == 3, "This function only exists for 3D elements");
  update_vertex_coord_(cell.vertex_coordinates());
  _qpoints = get_master_integration_points();
  _weights = get_master_integration_weights();
  update_();
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(const mesh::Cell & cell, const angem::Point<3,double> & point)
{
  static_assert(ELEMENT_DIM<vtk_id> == 3, "This function only exists for 3D elements");
  update_vertex_coord_(cell.vertex_coordinates());
  _weights = {1.0};
  _qpoints = {map_real_to_local_(point)};
  update_();
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(const mesh::Face & face, const angem::Point<3,double> & point)
{
  static_assert(ELEMENT_DIM<vtk_id> == 2, "This function only exists for 2D elements");
  update_vertex_coord_(face.vertex_coordinates());
  _weights = {1.0};
  _qpoints = {map_real_to_local_(point)};
  update_();
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(const mesh::Face & face)
{
  static_assert(ELEMENT_DIM<vtk_id> == 2, "This function only exists for 2D elements");
  update_vertex_coord_(face.vertex_coordinates());
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
        J( j, i ) += ref_grad[v][j] * _vertex_coord[v][i];
  return J;
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update_()
{
  _shape_values.resize( _qpoints.size() );
  _determinants.assign(_qpoints.size(), 0.0);
  _shape_grads.assign(_qpoints.size(), std::array<Point,N_ELEMENT_VERTICES<vtk_id>>());
  // update data in integration points
  for (std::size_t q=0; q<_qpoints.size(); ++q)
    update_(_qpoints[q], _shape_values[q], _shape_grads[q], _determinants[q]);
  // update values in center
  update_(_center, _shape_values_center, _shape_grads_center, _determinant_center);
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update_(const Point & p,
                               std::array<double,N_ELEMENT_VERTICES<vtk_id>> & shape_values,
                               std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & shape_grads,
                               double & determinant) const
{
  // compute shape values in q-points in master element
  update_shape_values_(p, shape_values);
  // compute shape grad in in master element
  const auto ref_grad = compute_ref_gradient_(p);
  angem::Tensor2<3, double> du_dx;
  if ( ELEMENT_DIM<vtk_id> == 3 )
    compute_detJ_and_invert_cell_jacobian_(ref_grad, du_dx, determinant);
  else  // dim == 2
    compute_detJ_and_invert_face_jacobian_(ref_grad, du_dx, determinant);
  // must be positive
  if ( determinant <= 0 ) throw std::runtime_error("Transformation Jacobian is not invertible");
  // compute the true shape function gradients
  update_shape_grads_(ref_grad, du_dx, shape_grads);
}


template<VTK_ID vtk_id>
void FeValues<vtk_id>::update_shape_grads_(const std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & ref_grad,
                                           const angem::Tensor2<3, double> & du_dx,
                                           std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & grads) const
{
    // compute the true shape function gradients
    for (size_t vertex = 0; vertex < N_ELEMENT_VERTICES<vtk_id>; ++vertex)
      for (size_t i=0; i<3; ++i)
        for (size_t j = 0; j < 3; ++j)
          // d phi_vert / dx_i = (d phi_vert / d u_j) * (d u_j / d x_i)
          grads[vertex][i] += ref_grad[vertex][j] * du_dx(j, i);
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update_vertex_coord_(const std::vector<Point> & vertex_coordinates)
{
  assert( vertex_coordinates.size() == N_ELEMENT_VERTICES<vtk_id> );
  for (size_t iv=0; iv<N_ELEMENT_VERTICES<vtk_id>; ++iv)
    _vertex_coord[iv] = vertex_coordinates[iv];
}

template<VTK_ID vtk_id>
double FeValues<vtk_id>::value(const size_t shape_index, const size_t qpoint) const
{
  assert( qpoint < _qpoints.size() && "qpoint too large" );
  assert( shape_index < N_ELEMENT_VERTICES<vtk_id> && "shape_index too large");
  return _shape_values[qpoint][shape_index];
}

template<VTK_ID vtk_id>
double FeValues<vtk_id>::value_center(const size_t shape_index) const
{
  assert( shape_index < N_ELEMENT_VERTICES<vtk_id> && "shape_index too large");
  return _shape_values_center[shape_index];
}

template<VTK_ID vtk_id>
Point FeValues<vtk_id>::grad(const size_t shape_index, const size_t qpoint) const
{
  assert( qpoint < _qpoints.size() && "qpoint too large" );
  assert( shape_index < N_ELEMENT_VERTICES<vtk_id> && "shape_index too large");
  return _shape_grads[qpoint][shape_index];
}

template<VTK_ID vtk_id>
Point FeValues<vtk_id>::grad_center(const size_t shape_index) const
{
  assert( shape_index < N_ELEMENT_VERTICES<vtk_id> && "shape_index too large");
  return _shape_grads_center[shape_index];
}

template<VTK_ID vtk_id>
double FeValues<vtk_id>::JxW(const size_t qpoint) const
{
  assert( qpoint < _qpoints.size() && "qpoint too large" );
  return _determinants[qpoint] * _weights[qpoint];
}

template<VTK_ID vtk_id>
double FeValues<vtk_id>::detJ_center() const
{
  return _determinant_center;
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::
update_shape_values_(const Point                                   & p,
                     std::array<double,N_ELEMENT_VERTICES<vtk_id>> & values) const
{
  for (size_t v=0; v<N_ELEMENT_VERTICES<vtk_id>; ++v)
    values[v] = eval_(p, v);
}

template<VTK_ID vtk_id>
std::array<Point,N_ELEMENT_VERTICES<vtk_id>> FeValues<vtk_id>::
compute_ref_gradient_(const Point & p) const
{
  std::array<Point,N_ELEMENT_VERTICES<vtk_id>> ref_grad;
  for (size_t v=0; v<N_ELEMENT_VERTICES<vtk_id>; ++v)
    ref_grad[v] = eval_derivative_(p, v);
  return ref_grad;
}

template<VTK_ID vtk_id>
void
FeValues<vtk_id>::
compute_detJ_and_invert_cell_jacobian_(const std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & ref_grad,
                                       angem::Tensor2<3, double> & du_dx,
                                       double & detJ) const
{
  // compute transformation jacobian
  const angem::Tensor2<3, double> dx_du = compute_jacobian_(ref_grad);
  // compute the determinant of transformation jacobian
  detJ = det(dx_du);
  // invert the jacobian to compute shape function gradients
  du_dx = invert(dx_du);
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::
compute_detJ_and_invert_face_jacobian_(const std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & ref_grad,
                                       angem::Tensor2<3, double> & J_inv,
                                       double & detJ) const
{
  // angem::Tensor2<3,double>
  angem::Plane<double> plane (_vertex_coord[0], _vertex_coord[1], _vertex_coord[2]);
  // const auto & basis = plane.get_basis();
  std::vector<Point> loc_coord(_vertex_coord.size());
  for (size_t i=0; i<_vertex_coord.size(); ++i)
    loc_coord[i] = plane.local_coordinates(_vertex_coord[i]);

  angem::Tensor2<3, double> J;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      for (std::size_t v=0; v<N_ELEMENT_VERTICES<vtk_id>; ++v)
        J( j, i ) += ref_grad[v][j] * loc_coord[v][i];

  // since the shape is 2d, we need to cast this down in 2D
  // since the third row and column of J are zero
  const angem::Tensor2<2, double> J2 = {J(0, 0), J(0, 1),
                                        J(1, 0), J(1, 1)};
  const angem::Tensor2<2, double> J2_inv = invert(J2);

  // cast it back to 3D since that's what the code uses to compute shape gradients
  J_inv = {J2_inv(0,0), J2_inv(0,1), 0.0,
           J2_inv(1,0), J2_inv(1,1), 0.0,
           0.0,         0.0,         0.0};
  detJ = det(J2);
}

template<> constexpr size_t N_ELEMENT_VERTICES<VTK_ID::TriangleID> = 3;
template<> constexpr size_t ELEMENT_DIM<VTK_ID::TriangleID> = 2;

template<> constexpr size_t N_ELEMENT_VERTICES<VTK_ID::TetrahedronID> = 4;
template<> constexpr size_t ELEMENT_DIM<VTK_ID::TetrahedronID> = 3;

template<> constexpr size_t N_ELEMENT_VERTICES<VTK_ID::HexahedronID> = 8;
template<> constexpr size_t ELEMENT_DIM<VTK_ID::HexahedronID> = 3;

}  // end namespace discretization
