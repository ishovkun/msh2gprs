#pragma once

#include "mesh/Mesh.hpp"
#include "angem/VTK_ID.hpp"
#include "angem/Tensor2.hpp"

namespace discretization {

using angem::VTK_ID;
using Point = angem::Point<3,double>;
// template<VTK_ID id> constexpr size_t N_ELEMENT_VERTICES;
// template<VTK_ID id> constexpr size_t ELEMENT_DIM;
template<VTK_ID id>
struct ElementTraits
{
  static const size_t n_vertices;
  static const size_t dim;
};

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
   * Update the needed internal quantities in the cell (or face) assigned as a vector of vertices.
   * Use this method before using the fem quantities in the new cell during the loop.
   */
  void update(std::vector<Point> & element_vertices);
  /**
   * Update the needed internal quantities in the cell in the gived integration points.
   * The point coordinate is the real (not reference) coordinates of the integration point.
   */
  void update(const mesh::Cell & cell, const std::vector<angem::Point<3,double>> & points);
  /**
   * Update the needed internal quantities in the face element.
   * Use this method before using the fem quantities in the new face during the loop.
   */
  void update(const mesh::Face & face);
  /**
   * Update the needed quantities in the face in the gived integration point.
   * The point coordinate is the real (not reference) coordinates of the integration point.
   */
  void update(const mesh::Face & face, const angem::Point<3,double> & point);
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

  void set_basis(const angem::Basis<3,double> & basis) { _face_basis = basis; _basis_set = true; }

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
               std::array<double,ElementTraits<vtk_id>::n_vertices> & shape_values,
               std::array<Point,ElementTraits<vtk_id>::n_vertices> & shape_grads,
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
  // compute shape values in q-points in master element
  void update_shape_values_(const Point                                   & p,
                            std::array<double,ElementTraits<vtk_id>::n_vertices> & values) const;
  // compute shape grads in q-points in master element
  void update_shape_grads_(const std::array<Point,ElementTraits<vtk_id>::n_vertices> & ref_grad,
                           const angem::Tensor2<3, double> & du_dx,
                           std::array<Point,ElementTraits<vtk_id>::n_vertices> & grads) const;
  // compute the gradients of shape fucntions in qpoint in master element
  std::array<Point,ElementTraits<vtk_id>::n_vertices> compute_ref_gradient_(const Point &p) const;

  /**
   * Map from real coordinates xyz to local coordinates in the master element uvw
   */
  Point map_real_to_local_(const Point & p) const;

  /**
   * Compute cell jacobian, its determinant, and invert it.
   */
  void compute_detJ_and_invert_cell_jacobian_(const std::array<Point,ElementTraits<vtk_id>::n_vertices> & ref_grad,
                                              angem::Tensor2<3, double> & du_dx,
                                              double & detJ) const;
  /**
   * Compute face jacobian, its determinant, and invert it.
   * Need a separate method for it since the jacobian will be 2x2 instead of 3x3.
   */
  void compute_detJ_and_invert_face_jacobian_(const std::array<Point,ElementTraits<vtk_id>::n_vertices> & ref_grad,
                                              angem::Tensor2<3, double> & du_dx,
                                              double & detJ) const;

  /*  ---------------------------------Variables ---------------------------------- */

  std::vector<Point> _qpoints;  // gauss point coordinates in master element
  Point              _center;   // center coordinates in master element
  std::vector<std::array<double,ElementTraits<vtk_id>::n_vertices>> _shape_values;
  std::vector<std::array<Point,ElementTraits<vtk_id>::n_vertices>> _shape_grads;
  std::array<double,ElementTraits<vtk_id>::n_vertices> _shape_values_center;
  std::array<Point,ElementTraits<vtk_id>::n_vertices>  _shape_grads_center;
  std::vector<double> _weights;
  std::vector<double> _determinants;
  double              _determinant_center;
  std::array<Point,ElementTraits<vtk_id>::n_vertices> _vertex_coord;
  angem::Basis<3,double> _face_basis;
  bool                   _basis_set = false;  // true after set_basis()
  std::array<size_t,ElementTraits<vtk_id>::n_vertices>    _vertex_ordering;
};

template<VTK_ID vtk_id>
FeValues<vtk_id>::FeValues()
{}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(const mesh::Cell & cell)
{
  static_assert(ElementTraits<vtk_id>::dim == 3, "This function only exists for 3D elements");
  update_vertex_coord_(cell.vertex_coordinates());
  _qpoints = get_master_integration_points();
  _weights = get_master_integration_weights();
  update_();
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(const mesh::Cell & cell, const std::vector<angem::Point<3,double>> & points)
{
  static_assert(ElementTraits<vtk_id>::dim == 3, "This function only exists for 3D elements");
  update_vertex_coord_(cell.vertex_coordinates());
  // weight = measure / n_points
  _weights = get_master_integration_weights();
  const double measure = std::accumulate(_weights.begin(), _weights.end(), 0.0);
  _weights.resize(points.size(), 0);
  for (size_t i = 0; i < _weights.size(); ++i)
    _weights[i] = measure / _weights.size();
  _qpoints.resize( points.size() );
  for (size_t q=0; q<points.size(); ++q)
    _qpoints[q] = map_real_to_local_(points[q]);
  update_();
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(const mesh::Face & face, const angem::Point<3,double> & point)
{
  static_assert(ElementTraits<vtk_id>::dim == 2, "This function only exists for 2D elements");
  update_vertex_coord_(face.vertex_coordinates());
  _weights = {1.0};
  _qpoints = {map_real_to_local_(point)};
  update_();
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(std::vector<Point> & element_vertices)
{
  update_vertex_coord_(element_vertices);
  _qpoints = get_master_integration_points();
  _weights = get_master_integration_weights();
  update_();
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(const mesh::Face & face)
{
  static_assert(ElementTraits<vtk_id>::dim == 2, "This function only exists for 2D elements");
  update_vertex_coord_(face.vertex_coordinates());
  _qpoints = get_master_integration_points();
  _weights = get_master_integration_weights();
  update_();
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update_()
{
  _shape_values.resize( _qpoints.size() );
  _determinants.assign(_qpoints.size(), 0.0);
  _shape_grads.resize(_qpoints.size());
  // update data in integration points
  for (std::size_t q=0; q<_qpoints.size(); ++q)
    update_(_qpoints[q], _shape_values[q], _shape_grads[q], _determinants[q]);
  // update values in center
  update_(_center, _shape_values_center, _shape_grads_center, _determinant_center);
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update_(const Point & p,
                               std::array<double,ElementTraits<vtk_id>::n_vertices> & shape_values,
                               std::array<Point,ElementTraits<vtk_id>::n_vertices> & shape_grads,
                               double & determinant) const
{
  // compute shape values in q-points in master element
  update_shape_values_(p, shape_values);
  // compute shape grad in in master element
  const auto ref_grad = compute_ref_gradient_(p);
  angem::Tensor2<3, double> du_dx;
  if ( ElementTraits<vtk_id>::dim == 3 )
    compute_detJ_and_invert_cell_jacobian_(ref_grad, du_dx, determinant);
  else  // dim == 2
    compute_detJ_and_invert_face_jacobian_(ref_grad, du_dx, determinant);
  // must be positive
  if ( determinant <= 0 ) throw std::runtime_error("Transformation det(J) is negative " + std::to_string(determinant));
  // compute the true shape function gradients
  update_shape_grads_(ref_grad, du_dx, shape_grads);
}


template<VTK_ID vtk_id>
void FeValues<vtk_id>::update_shape_grads_(const std::array<Point,ElementTraits<vtk_id>::n_vertices> & ref_grad,
                                           const angem::Tensor2<3, double> & du_dx,
                                           std::array<Point,ElementTraits<vtk_id>::n_vertices> & grads) const
{
  for (size_t vertex = 0; vertex < ElementTraits<vtk_id>::n_vertices; ++vertex)
    grads[vertex] = {0.0, 0.0, 0.0};

  // compute the true shape function gradients
  for (size_t vertex = 0; vertex < ElementTraits<vtk_id>::n_vertices; ++vertex)
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        // d phi_vert / dx_i = (d phi_vert / d u_j) * (d u_j / d x_i)
        grads[vertex][i] += ref_grad[vertex][j] * du_dx(j, i);
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update_vertex_coord_(const std::vector<Point> & vertex_coordinates)
{
  assert( vertex_coordinates.size() == ElementTraits<vtk_id>::n_vertices );
  for (size_t iv=0; iv<ElementTraits<vtk_id>::n_vertices; ++iv)
    _vertex_coord[iv] = vertex_coordinates[iv];

  std::iota(_vertex_ordering.begin(), _vertex_ordering.end(), 0);
  // check if vertex numbering is consistent with the current basis
  if ( ElementTraits<vtk_id>::dim == 2 )
  {
    if (_basis_set)
    {
      angem::Plane<double> plane (_vertex_coord[0], _vertex_coord[1], _vertex_coord[2]);
      if ( plane.normal() * _face_basis[2] < 0)
      {
        // reorder vertices
        std::reverse(_vertex_coord.begin()    + 1, _vertex_coord.end());
        std::reverse(_vertex_ordering.begin() + 1, _vertex_ordering.end());
      }
    }
  }
}

template<VTK_ID vtk_id>
double FeValues<vtk_id>::value(const size_t shape_index, const size_t qpoint) const
{
  assert( qpoint < _qpoints.size() && "qpoint too large" );
  assert( shape_index < ElementTraits<vtk_id>::n_vertices && "shape_index too large");
  // return _shape_values[qpoint][shape_index];
  return _shape_values[qpoint][_vertex_ordering[shape_index]];
}

template<VTK_ID vtk_id>
double FeValues<vtk_id>::value_center(const size_t shape_index) const
{
  assert( shape_index < ElementTraits<vtk_id>::n_vertices && "shape_index too large");
  return _shape_values_center[_vertex_ordering[shape_index]];
}

template<VTK_ID vtk_id>
Point FeValues<vtk_id>::grad(const size_t shape_index, const size_t qpoint) const
{
  assert( qpoint < _qpoints.size() && "qpoint too large" );
  assert( shape_index < ElementTraits<vtk_id>::n_vertices && "shape_index too large");
  return _shape_grads[qpoint][_vertex_ordering[shape_index]];
}

template<VTK_ID vtk_id>
Point FeValues<vtk_id>::grad_center(const size_t shape_index) const
{
  assert( shape_index < ElementTraits<vtk_id>::n_vertices && "shape_index too large");
  return _shape_grads_center[_vertex_ordering[shape_index]];
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
                     std::array<double,ElementTraits<vtk_id>::n_vertices> & values) const
{
  for (size_t v = 0; v < ElementTraits<vtk_id>::n_vertices; ++v)
    values[v] = eval_(p, v);
}

template<VTK_ID vtk_id>
std::array<Point,ElementTraits<vtk_id>::n_vertices> FeValues<vtk_id>::
compute_ref_gradient_(const Point & p) const
{
  std::array<Point,ElementTraits<vtk_id>::n_vertices> ref_grad;
  for (size_t v=0; v < ElementTraits<vtk_id>::n_vertices; ++v)
    ref_grad[v] = eval_derivative_(p, v);
  return ref_grad;
}

template<VTK_ID vtk_id>
void
FeValues<vtk_id>::
compute_detJ_and_invert_cell_jacobian_(const std::array<Point,ElementTraits<vtk_id>::n_vertices> & ref_grad,
                                       angem::Tensor2<3, double> & du_dx,
                                       double & detJ) const
{
  // compute transformation jacobian
  // see Becker E., Carey G., Oden J. Finite elements. An Introduction Volume 1 1981
  // Eq. 5.3.6
  angem::Tensor2<3, double> dx_du;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      for (std::size_t v=0; v<ElementTraits<vtk_id>::n_vertices; ++v)
        dx_du( i, j ) += ref_grad[v][j] * _vertex_coord[v][i];
  // compute the determinant of transformation jacobian
  detJ = det(dx_du);
  // invert the jacobian to compute shape function gradients
  du_dx = invert(dx_du);
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::
compute_detJ_and_invert_face_jacobian_(const std::array<Point,ElementTraits<vtk_id>::n_vertices> & ref_grad,
                                       angem::Tensor2<3, double> & J_inv,
                                       double & detJ) const
{
  angem::Plane<double> plane (_vertex_coord[0], _vertex_coord[1], _vertex_coord[2]);
  if (_basis_set)
  {
    plane.set_origin(_vertex_coord[0]);
    plane.set_basis(_face_basis);
  }
  // get vertex coordinates in 2d face basis
  std::array<Point,ElementTraits<vtk_id>::n_vertices> loc_coord;
  for (size_t i=0; i < ElementTraits<vtk_id>::n_vertices; ++i)
  {
    loc_coord[i] = plane.local_coordinates(_vertex_coord[i]);
    assert( std::fabs(loc_coord[i][2]) < 1e-12 );
  }

  // dx_j / du_i = \sum_k (d psi_k / du_j) * x_k_i
  angem::Tensor2<3, double> J;
  for (std::size_t i=0; i<2; ++i)
    for (std::size_t j=0; j<2; ++j)
      for (std::size_t v=0; v<ElementTraits<vtk_id>::n_vertices; ++v)
        J(i, j) += ref_grad[v][j] * loc_coord[v][i];
  J(2, 2) = 1;
  J_inv = invert(J);
  detJ = det(J);
}

template<VTK_ID vtk_id>
Point  FeValues<vtk_id>::map_real_to_local_(const Point & x) const
{
  /* This is a generic implementation of the mapping for non-linear shape functions.
   * Certain linear elements (e.g. triangles and tetras) use faster custom method.
   *
   * Algorithm:
   * Knowing the coordinate x in the real element and the
   * coordinates of the element vertices x_v, we want to
   * find the coordinate \xi in the master element by solving:
   * x = ⅀ᵥ xᵥ * φᵥ(ξ)
   * [Eq. 5.3.2 Becker E., Carey G., Oden J. Finite elements. An Introduction Volume 1 1981]
   *
   * We do this by linearizing this equation as follows:
   * residual R = ⅀ᵥ xᵥ φᵥ(ξ) - x
   * jacobian Jᵢⱼ = ∂Rᵢ / ∂ξⱼ = ⅀ᵥ xᵥᵢ * ∂φᵥ(ξ) / ∂ξⱼ
   * which is a 3x3 linear system.
   *
   * We wrap this system into a Newton-Rapson loop.
   */

  Point xi = {0,0,0};  // initial guess
  Point R;
  angem::Tensor2<3, double> J;
  double error = 1;
  for (size_t iter = 0; iter < 10; ++iter)
  {
    R.set_zero();
    J.set_zero();
    // assemble residual
    for (size_t v = 0; v < ElementTraits<vtk_id>::n_vertices; ++v)
      R += _vertex_coord[v] * eval_(xi, v);
    R -= x;
    // std::cout << "R.norm() = " << R.norm() << " xi = " << xi  << std::endl;

    if (R.norm() < 1e-8)  // converged
      return xi;

    // assemble Jacobian
    for (size_t v = 0; v < ElementTraits<vtk_id>::n_vertices; ++v)
    {
      for (size_t i = 0; i < ElementTraits<vtk_id>::dim; ++i)
      {
        const auto dphi_dxi = eval_derivative_(xi, v);
        for (size_t j = 0; j < ElementTraits<vtk_id>::dim; ++j)
          J(i, j) += _vertex_coord[v][i] * dphi_dxi[j];
      }
    }

    // solve system
    xi += invert(J) * (-R);
  }

  throw std::runtime_error("real-to-master mapping did not converge");
}

// template<> constexpr size_t N_ELEMENT_VERTICES<VTK_ID::TriangleID> = 3;
// template<> constexpr size_t ELEMENT_DIM<VTK_ID::TriangleID> = 2;

// template<> constexpr size_t N_ELEMENT_VERTICES<VTK_ID::TetrahedronID> = 4;
// template<> constexpr size_t ELEMENT_DIM<VTK_ID::TetrahedronID> = 3;

// template<> constexpr size_t N_ELEMENT_VERTICES<VTK_ID::HexahedronID> = 8;
// template<> constexpr size_t ELEMENT_DIM<VTK_ID::HexahedronID> = 3;

// template<> constexpr size_t N_ELEMENT_VERTICES<VTK_ID::QuadrangleID> = 4;
// template<> constexpr size_t ELEMENT_DIM<VTK_ID::QuadrangleID> = 2;

// template<> constexpr size_t N_ELEMENT_VERTICES<VTK_ID::WedgeID> = 6;
// template<> constexpr size_t ELEMENT_DIM<VTK_ID::WedgeID> = 3;

template<>
struct ElementTraits<VTK_ID::TriangleID>
{static const size_t n_vertices = 3; static const size_t dim = 2;};

template<>
struct ElementTraits<VTK_ID::QuadrangleID>
{static const size_t n_vertices = 4; static const size_t dim = 2;};

template<>
struct ElementTraits<VTK_ID::TetrahedronID>
{static const size_t n_vertices = 4; static const size_t dim = 3;};

template<>
struct ElementTraits<VTK_ID::HexahedronID>
{static const size_t n_vertices = 8; static const size_t dim = 3;};

template<>
struct ElementTraits<VTK_ID::WedgeID>
{static const size_t n_vertices = 6; static const size_t dim = 3;};

}  // end namespace discretization
