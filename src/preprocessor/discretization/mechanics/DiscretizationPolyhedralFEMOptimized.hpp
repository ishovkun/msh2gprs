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
  // void build() override;

 protected:
  void enumerate_elements_();
  void build_(mesh::Cell & cell) override;
  void reorder_faces_(mesh::Cell & dst, mesh::Cell const & src) const;
  // std::vector<uint64_t> compress_faces_(mesh::Cell const & cell) const;
  std::vector<uint64_t> compress_faces_(angem::Polyhedron<double> const & cell) const;
  std::vector<size_t> get_face_ordering_(angem::Polyhedron<double> const & dst,
                                         angem::Polyhedron<double> const & src) const;

  // returns element type; otherwise returns _shapes.size()
  size_t known_element_(angem::Polyhedron<double> const & poly, std::vector<size_t> & order) const;

  FiniteElementData scale_cell_fem_data_(mesh::Cell const & cell,
                                         FiniteElementData const & master_data);

  // *********** ATTRIBUTES ***********

  // a container for reference elements FEM data
  // std::unordered_map<size_t, std::vector<std::shared_ptr<PolyhedralElementBase>>> _masters;
  std::unordered_map<size_t, std::vector<size_t>> _vertex_to_cell_types;
  std::vector<std::shared_ptr<angem::Polyhedron<double>>> _shapes;
  std::vector<size_t> _element_types;
  std::vector<std::shared_ptr<PolyhedralElementBase>> _master_elements;
};

}  // end namespace discretization
