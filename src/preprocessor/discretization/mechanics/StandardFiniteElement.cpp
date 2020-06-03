#include "StandardFiniteElement.hpp"
#include "angem/VTK_ID.hpp"

namespace discretization {

using VTK_ID = angem::VTK_ID;

StandardFiniteElement::StandardFiniteElement(const mesh::Cell & cell,
                                             const bool update_face_values,
                                             const bool update_fracture_values)
    : _cell(cell), _update_frac_values(update_fracture_values)
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

  if (update_face_values || update_fracture_values)
    update_face_values_();

  if (update_fracture_values)
    update_cell_values_in_faces_();
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

void StandardFiniteElement::update_cell_values_in_faces_()
{
  const auto faces = _cell.faces();
  const size_t n_faces = faces.size();
  _face_fracture_data.resize( n_faces );

  const VTK_ID cell_id = static_cast<VTK_ID>(_cell.vtk_id());

  for (size_t i = 0; i < n_faces; ++i)
  {
    std::vector<Point> face_qpoints =
        compute_face_qpoint_coordinates(faces[i]->vertex_coordinates(), _face_data[i]);

    switch (cell_id)
    {
      case VTK_ID::TetrahedronID:
        {
          FeValues<VTK_ID::TetrahedronID> fe_values;
          fe_values.update(_cell, face_qpoints);
          build_<VTK_ID::TetrahedronID>(fe_values, _face_fracture_data[i], /*update_center = */ false);
          break;
        }
      case VTK_ID::HexahedronID:
        {
          FeValues<VTK_ID::HexahedronID> fe_values;
          fe_values.update(_cell, face_qpoints);
          build_<VTK_ID::HexahedronID>(fe_values, _face_fracture_data[i],
                                       /* update_center = */ false);
          break;
        }
      case VTK_ID::WedgeID:
        {
          FeValues<VTK_ID::WedgeID> fe_values;
          fe_values.update(_cell, face_qpoints);
          build_<VTK_ID::WedgeID>(fe_values, _face_fracture_data[i],
                                       /* update_center = */ false);
          break;
        }
      default:
        throw std::invalid_argument("FEM not implemented for this vtk type " + std::to_string(cell_id));
    }
  }

 
}

void StandardFiniteElement::update_face_values_()
{
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
