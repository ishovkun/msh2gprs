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
           const int                  marker,
           const std::size_t          parent )
    : m_index(cell_index),
      m_vertices(vertices),
      m_faces(faces),
      m_grid_vertices(grid_vertices),
      m_grid_cells(grid_cells),
      m_grid_faces(grid_faces),
      m_vtk_id(vtk_id),
      m_marker(marker),
      m_parent(parent)
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


std::vector<Cell*> Cell::neighbors()
{
  std::vector<size_t> neighbors_indices;
  neighbors_indices.reserve(m_faces.size());
  for (const size_t iface : m_faces)
    for ( const Cell * neighbor : m_grid_faces[iface].neighbors() )
      if (neighbor->index() != index())
        if (std::find(neighbors_indices.begin(), neighbors_indices.end(),
                      neighbor->index()) == neighbors_indices.end())
          neighbors_indices.push_back(neighbor->index());
  std::vector<Cell*> neighbor_cells;
  neighbor_cells.reserve(neighbors_indices.size());
  for (const size_t icell : neighbors_indices)
    neighbor_cells.push_back( &(m_grid_cells[icell]) );
  return neighbor_cells;
}


std::vector<const Cell*> Cell::neighbors() const
{
  std::vector<size_t> neighbors_indices;
  neighbors_indices.reserve(m_faces.size());
  for (const size_t iface : m_faces)
  {
    for ( const Cell * neighbor : m_grid_faces[iface].neighbors() )
    {
      if (neighbor->index() != index())
        if (std::find(neighbors_indices.begin(), neighbors_indices.end(),
                      neighbor->index()) == neighbors_indices.end())
          neighbors_indices.push_back(neighbor->index());
    }
  }
  std::vector<const Cell*> neighbor_cells;
  neighbor_cells.reserve(neighbors_indices.size());
  for (const size_t icell : neighbors_indices)
    neighbor_cells.push_back( &(m_grid_cells[icell]) );
  return neighbor_cells;
}



}  // end namespace
