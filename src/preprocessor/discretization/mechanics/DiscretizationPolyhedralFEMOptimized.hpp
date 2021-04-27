#pragma once
#include "DiscretizationPolyhedralFEM.hpp"
#include "Elements/PolyhedralElementBase.hpp"
#include "angem/Tensor2.hpp"
#include <memory>  // shared_ptr

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
  void reorder_faces_(mesh::Cell & dst, mesh::Cell const & src) const;
  std::vector<uint64_t> compress_faces_(mesh::Cell const & cell) const;

  // returns master if found; otherwise returns nullptr
  std::shared_ptr<PolyhedralElementBase>
  known_element_(mesh::Cell const & cell, std::vector<size_t> & order) const;

  FiniteElementData scale_cell_fem_data_(mesh::Cell const & cell,
                                         FiniteElementData const & master_data);

  // *********** ATTRIBUTES ***********

  // a container for reference elements FEM data
  std::unordered_map<size_t, std::vector<std::shared_ptr<PolyhedralElementBase>>> _masters;
};

}  // end namespace discretization
