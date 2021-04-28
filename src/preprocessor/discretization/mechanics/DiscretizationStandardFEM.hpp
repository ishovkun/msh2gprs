#pragma once
#include "DiscretizationFEMBase.hpp"

namespace discretization {

class DiscretizationStandardFEM : public DiscretizationFEMBase {
 public:
  DiscretizationStandardFEM(mesh::Mesh & grid,
                            const FiniteElementConfig & config,
                            const std::vector<int> & fracture_face_markers,
                            const std::vector<size_t> & neumann_face_indices);
  virtual ~DiscretizationStandardFEM() = default;

 protected:
  virtual void build_(mesh::Cell & cell) override;
};

}  // end namespace discretization
