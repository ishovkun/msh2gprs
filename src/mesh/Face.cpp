#include "Face.hpp"
#include "Cell.hpp"

namespace mesh
{

Face::Face(const std::size_t                       face_index,
           const std::size_t                       master_face_index,
           const std::vector<std::size_t> &        face_vertices,
           const int                               face_vtk_id,
           const int                               face_marker,
           std::vector<Cell>              &        grid_cells,
           std::vector<Point> &                    grid_vertices,
           std::vector<std::vector<std::size_t>> & grid_vertex_cells)
    : m_index(face_index),
      m_master_face_index(master_face_index),
      m_vertices(face_vertices),
      m_vtk_id(face_vtk_id),
      m_marker(face_marker),
      m_grid_cells(grid_cells),
      m_grid_vertices(grid_vertices),
      m_grid_vertex_cells(grid_vertex_cells)
{}


std::vector<Cell*> Face::neighbors()
{
  std::vector<size_t> face_neighbor_indices;
  face_neighbor_indices.reserve(2);
  for (const std::size_t vertex : vertices())
    for (const std::size_t cell_index : m_grid_vertex_cells[vertex])
    {
      auto & cell = m_grid_cells[cell_index];
      for (const Face * f : cell.faces())
        if (f->index() == index())
          if (std::find( face_neighbor_indices.begin(),
                         face_neighbor_indices.end(),
                         cell.index()) == face_neighbor_indices.end())
        {
          face_neighbor_indices.push_back(cell.index());
        }
    }
  std::vector<Cell*> face_neighbors;
  face_neighbors.reserve(face_neighbor_indices.size());
  for (const size_t icell : face_neighbor_indices)
    face_neighbors.push_back( &(m_grid_cells[icell]) );
  return face_neighbors;
}


std::vector<const Cell*> Face::neighbors() const
{
  std::vector<size_t> face_neighbor_indices;
  face_neighbor_indices.reserve(2);
  for (const std::size_t vertex_index : vertices())
    for (const std::size_t cell_index : m_grid_vertex_cells[vertex_index])
    {
      auto & cell = m_grid_cells[cell_index];
      for (const Face * f : cell.faces())
        if (f->index() == index())
          if (std::find( face_neighbor_indices.begin(),
                         face_neighbor_indices.end(),
                         cell.index()) == face_neighbor_indices.end())
        {
          face_neighbor_indices.push_back(cell.index());
        }
    }
  std::vector<const Cell*> face_neighbors;
  face_neighbors.reserve(face_neighbor_indices.size());
  for (const size_t icell : face_neighbor_indices)
    face_neighbors.push_back( &(m_grid_cells[icell]) );
  return face_neighbors;
}



std::vector<Point> Face::vertex_coordinates() const
{
  std::vector<Point> coordinates;
  coordinates.reserve(m_vertices.size());
  for (const std::size_t vertex_index : vertices())
    coordinates.push_back( m_grid_vertices[vertex_index] );
  return coordinates;
}


Polygon Face::polygon() const
{
  return angem::Polygon<double>(vertex_coordinates());
}


Point Face::normal() const
{
  return polygon().normal();
}


Point Face::center() const
{
  return polygon().center();
}

bool Face::has_vertex(const std::size_t vertex_index) const
{
  if (std::find(m_vertices.begin(),m_vertices.end(), vertex_index) ==
      m_vertices.end())
    return false;
  else return true;
}

}
