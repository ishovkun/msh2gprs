#pragma once
#include "DiscretizationPolyhedralFEM.hpp"
#include "angem/Tensor2.hpp"

namespace discretization {

class DiscretizationPolyhedralFEMOptimized : public DiscretizationPolyhedralFEM {
 public:
  // constructor
  DiscretizationPolyhedralFEMOptimized(mesh::Mesh & grid,
                                       const FiniteElementConfig & config,
                                       const std::vector<int> & fracture_face_markers,
                                       const std::vector<int> & neumann_face_markers);
  // destructor
  virtual ~DiscretizationPolyhedralFEMOptimized() = default;

 protected:
  void build_(mesh::Cell & cell) override;
  void build_cell_data_(mesh::Cell const & cell) override;

  bool known_element_(mesh::Cell const & cell,
                      std::vector<size_t> & order,
                      FiniteElementData const *& p_master) const;
  void scale_cell_fem_data_(mesh::Cell const & cell,
                            FiniteElementData const & data);
  void compute_detJ_and_invert_cell_jacobian_(const std::vector<angem::Point<3,double>> & ref_grad,
                                              angem::Tensor2<3, double> & du_dx,
                                              double & detJ,
                                              std::vector<angem::Point<3,double>> const & vertex_coord) const;
  void update_shape_grads_(std::vector<angem::Point<3,double>> const & ref_grads,
                           angem::Tensor2<3, double> const & du_dx,
                           std::vector<angem::Point<3,double>> &shape_grads) const;

  std::unordered_map<size_t, std::vector<FiniteElementData>> _cell_data_compressed;
};

}  // end namespace discretization
