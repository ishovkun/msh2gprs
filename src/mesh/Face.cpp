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
  // find how many times each vertex neighboring cell is encountered
  std::unordered_map<size_t, size_t> cell_time;
  for (const size_t vertex : m_vertices)
  {
    for (const size_t cell_index: m_grid_vertex_cells[vertex])
    {
      auto it = cell_time.find(cell_index);
      if (it == cell_time.end())
        cell_time.insert({cell_index, 1});
      else
        it->second++;
    }
  }

  // cells with the count == n_vertices are the face neighbors
  std::vector<Cell*> face_neighbors;
  face_neighbors.reserve(2);
  for (const auto it : cell_time)
  {
    if (it.second == m_vertices.size())
    {
      face_neighbors.push_back( &(m_grid_cells[it.first]) );
    }
  }
  return face_neighbors;
}


std::vector<const Cell*> Face::neighbors() const
{
  // find how many times each vertex neighboring cell is encountered
  std::unordered_map<size_t, size_t> cell_time;
  for (const size_t vertex : m_vertices)
  {
    for (const size_t cell_index: m_grid_vertex_cells[vertex])
    {
      auto it = cell_time.find(cell_index);
      if (it == cell_time.end())
        cell_time.insert({cell_index, 1});
      else
        it->second++;
    }
  }

  // cells with the count == n_vertices are the face neighbors
  std::vector<const Cell*> face_neighbors;
  face_neighbors.reserve(2);
  for (const auto it : cell_time)
  {
    if (it.second == m_vertices.size())
      face_neighbors.push_back( &(m_grid_cells[it.first]) );
  }
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
