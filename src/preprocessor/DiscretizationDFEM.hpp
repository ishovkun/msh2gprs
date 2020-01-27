#pragma once

#include "mesh/Mesh.hpp"

#ifdef WITH_GMSH
#include <gmsh.h>
#endif

namespace gprs_data
{

/* This class implements the Discrete Finite Element method (DFEM)
 * The Idea is to discretiza grid cells into simple shapes and compute
 * the shape functions with MSRSB method
 * (so that the simulator can use them as regular Finite Element
 * shape funcitons). */
class DiscretizationDFEM
{
 public:
  DiscretizationDFEM(const mesh::Mesh & grid);
  void build();

 protected:
  void build_(const angem::Polyhedron<double> & cell);
  void build_grid_(const angem::Polyhedron<double> & cell) const;
  void build_shape_functions_();
  double compute_element_size_(const angem::Polyhedron<double> & cell) const;

  const mesh::Mesh & _grid;
};

}  // end namepsace gprs_data
