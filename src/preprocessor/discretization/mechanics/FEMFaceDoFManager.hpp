#pragma once

#include "../flow/DoFNumbering.hpp"
#include "Mesh.hpp"
#include <vector>

namespace discretization {

/* This class is a utility to number degrees of freedom for
 * a face linear system used in polyhedral FEM */
class FEMFaceDoFManager
{
 public:
  FEMFaceDoFManager() {};
  DoFNumbering build(const mesh::Mesh & grid,
                     const std::vector<size_t> & face_indices) const;
};

}  // end namespace discretization
