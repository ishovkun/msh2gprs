#pragma once
#include "DiscretizationPolyhedralFEM.hpp"

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
  bool known_element_(mesh::Cell const & cell,
                      std::vector<size_t> & order,
                      FiniteElementDataTopology const *& p_master) const;
  void scale_cell_fem_data_(mesh::Cell const & cell,
                            FiniteElementDataTopology const & data);

  std::unordered_map<size_t, std::vector<FiniteElementDataTopology>> _cell_data_compressed;
};

}  // end namespace discretization
