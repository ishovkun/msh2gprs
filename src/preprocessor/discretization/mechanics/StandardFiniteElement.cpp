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
      {
        FeValues<VTK_ID::TetrahedronID> fe_values;
        fe_values.update(_cell);
        build_<VTK_ID::TetrahedronID>(fe_values, _cell_data);
        break;
      }
    case VTK_ID::HexahedronID:
      {
        FeValues<VTK_ID::HexahedronID> fe_values;
        fe_values.update(_cell);
        build_<VTK_ID::HexahedronID>(fe_values, _cell_data);
        break;
      }
    default:
      throw std::invalid_argument("FEM not implemented for this vtk type");
  }

  const auto faces = _cell.faces();
  _face_data.resize( faces.size() );
  for (size_t i=0; i<faces.size(); ++i)
  {
    const mesh::Face* face = faces[i];
    const VTK_ID id = static_cast<VTK_ID>(face->vtk_id());
    switch (id)
    {
      case VTK_ID::TriangleID:
        {
          FeValues<VTK_ID::TriangleID> fe_values;
          fe_values.update(*face);
          build_<VTK_ID::TriangleID>(fe_values, _face_data[i]);
          break;
        }
      case VTK_ID::QuadrangleID:
        {
          FeValues<VTK_ID::QuadrangleID> fe_values;
          fe_values.update(*face);
          build_<VTK_ID::QuadrangleID>(fe_values, _face_data[i]);
          break;
        }
      default:
        throw std::invalid_argument("FEM not implemented for this vtk type");
    }
  }
}

}  // end namespace discretization
