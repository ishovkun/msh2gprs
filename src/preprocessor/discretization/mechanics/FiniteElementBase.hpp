#pragma once
#include "FiniteElementData.hpp"            // provides FiniteElementData

namespace discretization {

class FiniteElementBase
{
 public:
  // get FE data for volume integration
  const FiniteElementData & get_cell_data() const {return _cell_data;}
  // get FE data for surface integration
  const std::vector<FiniteElementData> & get_face_data() const {return _face_data;}

 protected:
  FiniteElementData _cell_data;                // FEM values and gradients in cell integration points
  std::vector<FiniteElementData> _face_data;   // FEM values and gradients in face integration points
};

}  // end namespace discretization
