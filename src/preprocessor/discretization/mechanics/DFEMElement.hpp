#pragma once

#include "mesh/Cell.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include <Eigen/Sparse>

namespace discretization {

class DFEMElement
{
 public:
  DFEMElement(const mesh::Cell & cell);
 protected:
  void build_();
  void build_shape_functions_();
  void build_jacobian_();
  // // fill out element_numbering and node_numbering
  void numberNodesEndElements_(std::vector<int> &element_types,
                        std::vector<std::vector<std::size_t> > & element_tags,
                        const std::vector<std::vector<std::size_t> > &node_tags);


 private:
  const mesh::Cell & _cell;
  std::unordered_map<size_t, size_t> _cell_numbering;
  std::unordered_map<size_t, size_t> _node_numbering;
};

}  // end namespace discretization
