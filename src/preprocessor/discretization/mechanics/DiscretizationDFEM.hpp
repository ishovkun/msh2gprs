#pragma once

#include "mesh/Mesh.hpp"
#include "FiniteElementData.hpp"

namespace discretization
{

/* This class implements the Discrete Finite Element method (DFEM)
 * The Idea is to discretiza grid cells into simple shapes and compute
 * the shape functions with MSRSB method
 * (so that the simulator can use them as regular Finite Element
 * shape funcitons). */
class DiscretizationDFEM
{
 public:
  DiscretizationDFEM(const mesh::Mesh & grid,
                     const double       msrsb_tol);
  void build();

  // get vector of finite element data that corresponds to 3D cells
  std::vector<FiniteElementData> get_cell_data() const;
  // get vector of finite element data that corresponds to faces of 3D cells
  std::vector<FiniteElementData> get_face_data() const;

 protected:
  const mesh::Mesh & _grid;
  const double _msrsb_tol;  // msrsb tolerance
};

}  // end namepsace discretization
