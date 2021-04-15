#pragma once
#include "DiscretizationFEMBase.hpp"

namespace discretization {

class DiscretizationPolyhedralFEM : public DiscretizationFEMBase {
 public:
  DiscretizationPolyhedralFEM(mesh::Mesh & grid,
                              const FiniteElementConfig & config,
                              const std::vector<int> & fracture_face_markers,
                              const std::vector<int> & neumann_face_markers);
  virtual ~DiscretizationPolyhedralFEM() = default;
 protected:
  void build_(mesh::Cell & cell) override;
};

}  // end namespace discretization
