#pragma once

#include "mesh/Cell.hpp"                    // provices mesh::Cell
#include "FiniteElementData.hpp"            // provides FiniteElementData

namespace discretization {

/**
 * This class implements a standard finite element discretization for a
 * single cell of the mesh.
 */
class StandardFiniteElement {
 public:
  /**
   * Constructor.
   * Build FEM discretization of the cell.
   * Input:
   * @param[in] cell : grid cell to be discretized
   */
  StandardFiniteElement(const mesh::Cell & cell);
  // get FE data for volume integration
  const FiniteElementData & get_cell_data() const {return _cell_data;}
  // get FE data for surface integration
  const std::vector<FiniteElementData> & get_face_data() const {return _face_data;}

 protected:
  void build_();

  const mesh::Cell & _cell;
  FiniteElementData _cell_data;                // FEM values and gradients in cell integration points
  std::vector<FiniteElementData> _face_data;   // FEM values and gradients in face integration points
};

}  // end namespace discretization
