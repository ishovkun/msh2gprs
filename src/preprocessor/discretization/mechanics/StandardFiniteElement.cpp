#include "StandardFiniteElement.hpp"
#include "angem/VTK_ID.hpp"

namespace discretization {

using VTK_ID = angem::VTK_ID;

StandardFiniteElement::StandardFiniteElement(const mesh::Cell & cell)
    : _cell(cell)
{
  const VTK_ID id = static_cast<VTK_ID>(_cell.vtk_id());
  switch (id)
  {
    case VTK_ID::TetrahedronID:
      build_<VTK_ID::TetrahedronID>();
      break;
    case VTK_ID::HexahedronID:
      build_<VTK_ID::HexahedronID>();
      break;
  }
}

}  // end namespace discretization
