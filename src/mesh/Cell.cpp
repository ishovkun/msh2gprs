#include "Cell.hpp"
#include "angem/PolyhedronFactory.hpp"

namespace mesh
{

Cell::Cell(const std::size_t          cell_index,
           const std::vector<std::size_t> & vertices,
           const std::vector<std::size_t> & faces,
           std::vector<Point>        & grid_vertices,
           std::vector<Cell>        & grid_cells,
           std::vector<Face>        & grid_faces,
           const int                  vtk_id,
           const int                  marker)
    : m_index(cell_index),
      m_vertices(vertices),
      m_faces(faces),
      m_grid_vertices(grid_vertices),
      m_grid_cells(grid_cells),
      m_grid_faces(grid_faces),
      m_vtk_id(vtk_id),
      m_marker(marker)
{}

std::unique_ptr<Polyhedron> Cell::polyhedron() const
{
  return angem::PolyhedronFactory::create<double>(m_grid_vertices,
                                                  m_vertices, vtk_id());
}


std::vector<Face*> Cell::faces()
{
  std::vector<Face*> cell_faces;   cell_faces.reserve(m_faces.size());
  for (const std::size_t face_index: m_faces)
    cell_faces.push_back(&m_grid_faces[face_index]);
  return cell_faces;
 }


std::vector<const Face*> Cell::faces() const
{
  std::vector<const Face*> cell_faces;
  cell_faces.reserve(m_faces.size());
  for (const std::size_t face_index: m_faces)
    cell_faces.push_back(&m_grid_faces[face_index]);
  return cell_faces;
}


double Cell::volume() const
{
  return polyhedron()->volume();
}


Point Cell::center() const
{
  return polyhedron()->center();
}


}  // end namespace
