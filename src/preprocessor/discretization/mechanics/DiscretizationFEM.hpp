#pragma once

#include "mesh/Mesh.hpp"
#include "FiniteElementData.hpp"
#include "FiniteElementBase.hpp"
#include "config/FiniteElementConfig.hpp"
#include "angem/Basis.hpp"

namespace discretization
{

struct FaceOrientation
{
  angem::Basis<3, double> basis;
  bool assigned;
};

/* This class implements Finite Element discretization
 * The Idea is to compute shape functions and feed them to the simulator
 */
class DiscretizationFEM
{
 public:
  DiscretizationFEM(const mesh::Mesh & grid, const FiniteElementConfig & config,
                    const std::vector<int> & fracture_face_markers,
                    const std::vector<int> & neumann_face_markers);

  // get vector of finite element data that corresponds to 3D cells
  const std::vector<FiniteElementData> & get_cell_data() const { return _cell_data; }
  // get vector of finite element data that corresponds to faces of 3D cells
  const std::vector<FiniteElementData> & get_face_data() const { return _face_data; }
  // get vector fe-fracture data (vector[face][cell])
  const std::vector<std::vector<FiniteElementData>> & get_fracture_data() const { return _frac_data; }

 protected:
  // just a debug function
  void analyze_cell_(const mesh::Cell & cell);
  // choose an element discretization based on config and cell vtk id
  std::unique_ptr<FiniteElementBase> build_element(const mesh::Cell & cell);

  const angem::Basis<3, double>& get_basis_(const mesh::Face & face,
                                            FaceOrientation &orientation) noexcept;

  const mesh::Mesh & _grid;
  const FiniteElementConfig & _config;
  std::unordered_map<int, FaceOrientation> _fracture_face_orientation;
  std::unordered_map<int, FaceOrientation> _neumann_face_orientation;
  // std::unordered_set<int> _fracture_face_markers;
  // std::unordered_set<int> _neumann_face_markers;
  std::vector<FiniteElementData> _cell_data, _face_data;
  std::vector<std::vector<FiniteElementData>> _frac_data;
};

}  // end namepsace discretization
