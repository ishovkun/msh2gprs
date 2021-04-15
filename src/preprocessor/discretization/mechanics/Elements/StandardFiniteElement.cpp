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
    case VTK_ID::WedgeID:
      {
        FeValues<VTK_ID::WedgeID> fe_values;
        fe_values.update(_cell);
        build_<VTK_ID::WedgeID>(fe_values, _cell_data);
        break;
      }
    default:
      throw std::invalid_argument("FEM not implemented for this vtk type " + std::to_string(id));
  }
}

std::vector<angem::Point<3,double>> compute_face_qpoint_coordinates(const std::vector<angem::Point<3,double>> & face_vertices,
                                                                    const FiniteElementData & fe_data)
{
  std::vector<Point> face_qpoints(fe_data.points.size());

  for (size_t v = 0; v < face_vertices.size(); ++v)
    for (size_t q = 0; q < face_qpoints.size(); ++q)
      face_qpoints[q] += face_vertices[v] * fe_data.points[q].values[v];

  return face_qpoints;
}

FiniteElementData StandardFiniteElement::get_fracture_data(const size_t iface,
                                                           const angem::Basis<3,double> basis)
{
  const auto faces = _cell.faces();
  const size_t n_faces = faces.size();

  assert( iface < n_faces );

  FiniteElementData data;   // FEM values and gradients in face integration points

  const auto face_qpoints = compute_face_qpoint_coordinates(faces[iface]->vertex_coordinates(),
                                                            get_face_data(iface, basis));

  const VTK_ID cell_id = static_cast<VTK_ID>(_cell.vtk_id());
    switch (cell_id)
    {
      case VTK_ID::TetrahedronID:
        {
          FeValues<VTK_ID::TetrahedronID> fe_values;
          fe_values.update(_cell, face_qpoints);
          build_<VTK_ID::TetrahedronID>(fe_values, data, /*update_center = */ false);
          break;
        }
      case VTK_ID::HexahedronID:
        {
          FeValues<VTK_ID::HexahedronID> fe_values;
          fe_values.update(_cell, face_qpoints);
          build_<VTK_ID::HexahedronID>(fe_values, data, /* update_center = */ false);
          break;
        }
      case VTK_ID::WedgeID:
        {
          FeValues<VTK_ID::WedgeID> fe_values;
          fe_values.update(_cell, face_qpoints);
          build_<VTK_ID::WedgeID>(fe_values, data, /* update_center = */ false);
          break;
        }
      default:
        throw std::invalid_argument("FEM not implemented for this vtk type " + std::to_string(cell_id));
    }

    return data;
}

FiniteElementData StandardFiniteElement::get_face_data(const size_t iface,
                                                       const angem::Basis<3,double> basis)
{
  const auto faces = _cell.faces();
  FiniteElementData data;

  const mesh::Face* face = faces[iface];

  const VTK_ID id = static_cast<VTK_ID>(face->vtk_id());
  switch (id)
  {
    case VTK_ID::TriangleID:
      {
        FeValues<VTK_ID::TriangleID> fe_values;
        fe_values.update(*face);
        fe_values.set_basis(basis);
        build_<VTK_ID::TriangleID>(fe_values, data);
        break;
      }
    case VTK_ID::QuadrangleID:
      {
        FeValues<VTK_ID::QuadrangleID> fe_values;
        fe_values.update(*face);
        fe_values.set_basis(basis);
        build_<VTK_ID::QuadrangleID>(fe_values, data);
        break;
      }
    default:
      throw std::invalid_argument("FEM not implemented for this vtk type");
  }
  return data;
}


}  // end namespace discretization
