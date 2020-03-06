#pragma once
#include "angem/Point.hpp"
#include <vector>

namespace gprs_data {

using Point = angem::Point<3,double>;

class FeValues
{
 public:
  // Constructor
  FeValues(const int element_type, const size_t n_elements);
  /* update function gradients and jacobians with
   * custom points inside the element.
   * Note that the points must be inside the element*/
  void update(const size_t element_tag, const std::vector<Point> & points );
  /* update function gradients and jacobians */
  void update(const size_t element_tag);
  // get value of shape i in integration point q
  double value(const size_t i, const size_t q) const;
  // get gradient of shape i in integration point q
  Point grad(const size_t i, const size_t q) const;
  // get determinant
  double JxW(const size_t q) const;
  // number of integration points
  size_t n_q_points() const;
  // number of vertices in the reference element
  size_t n_vertices() const;

 protected:
  void initialize_();
  void compute_shape_grads_();
  void get_elements_();

 private:
  void debug_print_cell_config();
  const int _element_type;
  const size_t _n_elements;
  int _n_comp;
  size_t _cell_index;
  std::vector<size_t> _element_tags;
  std::vector<double> _ref_points;   // integration on reference element
  std::vector<double> _weights;      // integration weights
  // basis function gradients [dxi_duj] on ref element
  std::vector<double> _ref_values;
  // derivatives in the reference element: d (shape) / du
  std::vector<double> _ref_gradients;
  // jacobians
  // [e1g1Jxu, e1g1Jyu, e1g1Jzu, e1g1Jxv, ..., e1g1Jzw, e1g2Jxu, ..., e1gGJzw, e2g1Jxu, ...]
  std::vector<double> _all_jacobians;
  // determinant of the Jacobian matrix at each integration
  // point: [e1g1, e1g2, ... e1gG, e2g1, ...]
  std::vector<double> _all_determinants;
  std::vector<double> _true_points;  // integration points on real element
  // determinant values in each gauss point for the current cell
  std::vector<double> _determinants;
  // for only a single element at a time
  std::vector<double> _grad;  // shape gradients on real element
  std::vector<double> _inv_determinants;
};

}  // end namespace gprs_data
