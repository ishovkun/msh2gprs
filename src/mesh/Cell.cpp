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
      pm_grid_vertices(&grid_vertices),
      pm_grid_cells(&grid_cells),
      pm_grid_faces(&grid_faces),
      m_vtk_id(vtk_id),
      m_marker(marker),
      m_parent(parent)
{
  if (m_parent == constants::invalid_index)
    m_parent = m_index;
}

std::unique_ptr<Polyhedron> Cell::polyhedron() const
{
  std::vector<std::vector<std::size_t>> faces_vertices;
  for (const Face* face : faces())
    faces_vertices.push_back(face->vertices());
  return std::make_unique<angem::Polyhedron<double>>(*pm_grid_vertices,
                                                     faces_vertices, vtk_id());
}


std::vector<Face*> Cell::faces()
{
  std::vector<Face*> cell_faces;   cell_faces.reserve(m_faces.size());
  for (const std::size_t face_index: m_faces)
    cell_faces.push_back(&(*pm_grid_faces)[face_index]);
  return cell_faces;
 }


std::vector<const Face*> Cell::faces() const
{
  std::vector<const Face*> cell_faces;
  cell_faces.reserve(m_faces.size());
  for (const std::size_t face_index: m_faces)
    cell_faces.push_back(&(*pm_grid_faces)[face_index]);
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
    for ( const Cell * neighbor : (*pm_grid_faces)[iface].neighbors() )
      if (neighbor->index() != index())
        if (std::find(neighbors_indices.begin(), neighbors_indices.end(),
                      neighbor->index()) == neighbors_indices.end())
          neighbors_indices.push_back(neighbor->index());
  std::vector<Cell*> neighbor_cells;
  neighbor_cells.reserve(neighbors_indices.size());
  for (const size_t icell : neighbors_indices)
    neighbor_cells.push_back( &((*pm_grid_cells)[icell]) );
  return neighbor_cells;
}


std::vector<const Cell*> Cell::neighbors() const
{
  std::vector<size_t> neighbors_indices;
  neighbors_indices.reserve(m_faces.size());
  for (const size_t iface : m_faces)
  {
    for ( const Cell * neighbor : (*pm_grid_faces)[iface].neighbors() )
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
    neighbor_cells.push_back( &((*pm_grid_cells)[icell]) );
  return neighbor_cells;
}

Cell & Cell::operator=(const Cell & other)
{
  m_index = other.index();
  m_vertices = other.vertices();
  m_faces = other.m_faces;
  pm_grid_vertices = other.pm_grid_vertices;
  pm_grid_cells = other.pm_grid_cells;
  pm_grid_faces = other.pm_grid_faces;
  m_vtk_id = other.vtk_id();
  m_marker = other.marker();
  m_children = other.children();
  return *this;
}

const Cell & Cell::ultimate_parent() const
{
  const Cell * par = &parent();
  while (par->parent().index() != par->index())
    par = &par->parent();
  return *par;
}

Cell & Cell::ultimate_parent()
{
  return const_cast<Cell&>(ultimate_parent());
}


std::vector<size_t> Cell::ultimate_children() const
{
  assert( false && "write proper code for ultimate children" );
  std::vector<size_t>  ch = children();
  for (const size_t child_index : ch)
  {
    const Cell & child = (*pm_grid_cells)[child_index];
    for (const size_t grand_child : child.ultimate_children())
      ch.push_back( grand_child );
  }
  return ch;
}

}  // end namespace
