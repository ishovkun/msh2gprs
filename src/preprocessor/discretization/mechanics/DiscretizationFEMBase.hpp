#pragma once

#include "mesh/Mesh.hpp"
#include "mesh/Cell.hpp"
#include "Elements/FiniteElementBase.hpp"
#include "DiscretizationFEMBase.hpp"
#include "config/FiniteElementConfig.hpp"
#include "angem/Basis.hpp"


namespace discretization {

struct FaceOrientation
{
  angem::Basis<3, double> basis;
  bool assigned;
};

class DiscretizationFEMBase {
 public:
  virtual void build();

  // get vector of finite element data that corresponds to 3D cells
  const std::vector<FiniteElementData> & get_cell_data() const { return _cell_data; }
  // get vector of finite element data that corresponds to faces of 3D cells
  const std::vector<FiniteElementData> & get_face_data() const { return _face_data; }
  // get vector fe-fracture data (vector[face][cell])
  const std::vector<std::vector<FiniteElementData>> & get_fracture_data() const { return _frac_data; }
  // destructor
  virtual ~DiscretizationFEMBase() = default;

 protected:
  DiscretizationFEMBase(mesh::Mesh & grid, const FiniteElementConfig & config,
                        const std::vector<int> & fracture_face_markers,
                        const std::vector<int> & neumann_face_markers);
  virtual void build_(mesh::Cell & cell) = 0;
  virtual void build_cell_data_(mesh::Cell const & cell);
  void build_face_data_(mesh::Cell const & cell);
  const angem::Basis<3, double>& get_basis_(const mesh::Face & face,
                                            FaceOrientation &orientation) noexcept;

 protected:
  mesh::Mesh & _grid;  // non-constant since there could be local vertex reordering
  const FiniteElementConfig & _config;
  std::unordered_map<int, FaceOrientation> _fracture_face_orientation;
  std::unordered_map<int, FaceOrientation> _neumann_face_orientation;
  std::vector<FiniteElementData> _cell_data, _face_data;
  std::vector<std::vector<FiniteElementData>> _frac_data;
  std::unordered_map<size_t, std::vector<FiniteElementData>> _cell_data_compressed;
  angem::Basis<3, double> _face_basis;
  std::unique_ptr<FiniteElementBase> _element;
};

}  // end namespace discretization
