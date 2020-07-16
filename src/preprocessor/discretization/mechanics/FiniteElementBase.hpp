#pragma once
#include "FiniteElementData.hpp"            // provides FiniteElementData
#include "angem/Basis.hpp"

namespace discretization {

class FiniteElementBase
{
 public:
  // get FE data for volume integration
  virtual FiniteElementData get_cell_data() {return _cell_data;}
  // get FE data for surface integration
  virtual FiniteElementData get_face_data(const size_t iface,
                                          const angem::Basis<3,double> basis) = 0;
  // get FE data of cell shape functions at face integration points.
  // This is needed for modeling discrete fractures
  virtual FiniteElementData get_fracture_data(const size_t iface,
                                              const angem::Basis<3,double> basis) = 0;
  // default destructor
  virtual ~FiniteElementBase() = default;

 protected:
  FiniteElementData _cell_data;                // FEM values and gradients in cell integration points
};

}  // end namespace discretization
