#pragma once
#include "FiniteElementData.hpp"            // provides FiniteElementData
#include "angem/Basis.hpp"
#include "mesh/Cell.hpp"
#include "mesh/Face.hpp"

namespace discretization {

class FiniteElementBase
{
 public:
  // get FE data for volume integration
  virtual FiniteElementData get_cell_data() {return _cell_data;}
  // get FE data for surface integration
  virtual FiniteElementData get_face_data(const size_t iface) = 0;
  // get FE data of cell shape functions at face integration points.
  // This is needed for modeling discrete fractures
  virtual FiniteElementData get_fracture_data(const size_t iface,
                                              const angem::Basis<3,double> basis) = 0;
  // default destructor
  virtual ~FiniteElementBase() = default;

 protected:
  // Get the right-hand basis of the face (that belongs to cell)
  // so that the normal vector (component 3) is orented out of the
  // cell
  static angem::Basis<3,double> get_face_basis_(mesh::Face const & face,
                                                mesh::Cell const & cell);

  FiniteElementData _cell_data;                // FEM values and gradients in cell integration points
};

}  // end namespace discretization
