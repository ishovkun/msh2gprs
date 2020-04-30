#pragma once

#include "mesh/Mesh.hpp"
#include "FiniteElementData.hpp"
#include "FiniteElementBase.hpp"
#include "config/FiniteElementConfig.hpp"

namespace discretization
{

/* This class implements Finite Element discretization
 * The Idea is to compute shape functions and feed them to the simulator
 */
class DiscretizationFEM
{
 public:
  DiscretizationFEM(const mesh::Mesh & grid, const FiniteElementConfig & config,
                    const std::vector<int> & fracture_face_markers);

  // get vector of finite element data that corresponds to 3D cells
  const std::vector<FiniteElementData> & get_cell_data() const { return _cell_data; }
  // get vector of finite element data that corresponds to faces of 3D cells
  const std::vector<FiniteElementData> & get_face_data() const { return _face_data; }

 protected:
  // just a debug function
  void analyze_cell_(const mesh::Cell & cell);
  // choose an element discretization based on config and cell vtk id
  std::unique_ptr<FiniteElementBase> build_element(const mesh::Cell & cell);
  // true if at least one of the cell faces is a fracture
  bool has_fracture_face_(const mesh::Cell & cell);

  const mesh::Mesh & _grid;
  const FiniteElementConfig & _config;
  std::vector<FiniteElementData> _cell_data, _face_data;
  std::unordered_set<int> _fracture_face_markers;
};

}  // end namepsace discretization
