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
                                       const std::vector<size_t> & neumann_face_indices);
  // destructor
  virtual ~DiscretizationPolyhedralFEMOptimized() = default;

  std::vector<int> get_cell_isomorphic_groups() const override;

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
  std::unordered_map<size_t, std::vector<std::shared_ptr<angem::Polyhedron<double>>>> _shapes;
  std::vector<int> _groups;
  size_t _ngroups{0};
};

}  // end namespace discretization
