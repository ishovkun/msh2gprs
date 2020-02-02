#pragma once

#include "mesh/Mesh.hpp"
#include "gmsh_interface/GmshInterface.hpp"

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
  // build discretization for a single cell of the original grid
  // that consists of meshing the element and computing shape functons
  void build_(const mesh::Cell & cell);
  // build dfem shape functions
  void build_shape_functions_();
  // build the jacobian of a laplace equation of a single element of the original
  // mesh that is discretized with gmsh
  void build_jacobian_();
  void build_local_matrix_(const int element_type, const size_t element_tag);
  // fill out element_numbering and node_numbering
  void numberNodesEndElements_(std::vector<int> &element_types,
                        std::vector<std::vector<std::size_t> > & element_tags,
                        const std::vector<std::vector<std::size_t> > &node_tags);

  std::unordered_map<size_t, size_t> _cell_numbering;
  std::unordered_map<size_t, size_t> _node_numbering;
  const mesh::Mesh & _grid;
};

}  // end namepsace gprs_data
