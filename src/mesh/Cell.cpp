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
  std::vector<angem::Point<3,double>> coord = vertex_coordinates();
  const auto verts = vertices();
  std::map<size_t, size_t> old_to_new;
  size_t i = 0;
  for (auto v : verts)
    old_to_new[v] = i++;

  std::vector<std::vector<size_t>> poly_faces;
  for (auto face : faces())
  {
    std::vector<size_t> face_vertices;
    for (auto v : face->vertices())
    {
      const size_t idx = old_to_new[v];
      face_vertices.push_back(idx);
    }
    poly_faces.push_back(std::move(face_vertices));
  }
  const angem::Polyhedron<double> poly(coord, poly_faces);
  return std::make_unique<angem::Polyhedron<double>>(coord, poly_faces, vtk_id());
}

std::vector<Point> Cell::vertex_coordinates() const
{
  std::vector<Point> coordinates;
  coordinates.reserve(m_vertices.size());
  for (const std::size_t vertex_index : vertices())
    coordinates.push_back( (*pm_grid_vertices)[vertex_index] );
  return coordinates;
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
  m_children = other.m_children;
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
  const auto & const_this = *this;
  return const_cast<Cell&>(const_this.ultimate_parent());
}


std::vector<const Cell*> Cell::immediate_children() const
{
  std::vector<const Cell*> result;
  for (const size_t ichild : m_children)
    result.push_back( &((*pm_grid_cells)[ichild]) );
  return result;
}


std::vector<const Cell*> Cell::ultimate_children() const
{
  std::vector<size_t>  uc;
  this->ultimate_children_(uc);
  std::vector<const Cell*> result;
  for (const size_t ichild : uc)
    result.push_back( &((*pm_grid_cells)[ichild]) );
  return result;
}

void Cell::ultimate_children_(std::vector<size_t> & uc) const
{
  if ( m_children.empty() )
    uc.push_back(index());
  else
  {
    for (const Cell * p_child : immediate_children())
      p_child->ultimate_children_(uc);
  }
}

std::vector<const Cell *> Cell::all_level_children() const
{
  std::vector<size_t> ichildren;
  all_level_children_(ichildren);
  std::vector<const Cell*> result;
  for (const size_t ichild : ichildren)
    result.push_back( &((*pm_grid_cells)[ichild]) );
  return result;
}

void Cell::all_level_children_(std::vector<size_t> & ichildren) const
{
  ichildren.push_back(index());
  for (const Cell * p_child : immediate_children())
    p_child->all_level_children_(ichildren);
}

bool Cell::has_edge(const vertex_pair edge) const
{
  for (auto face : faces())
    if (face->has_edge(edge))
      return true;
  return false;
}

bool Cell::has_vertex(const size_t vert) const
{
  for (const auto & face : faces())
    if (face->has_vertex(vert))
      return true;
  return false;
}

// std::vector<size_t> Cell::sorted_vertices() const
// {
//   std::vector<size_t> result(vertices());
//   std::sort(result.begin(), result.end());
//   return result;
// }

std::vector<vertex_pair> Cell::edges() const noexcept
{
  std::set<vertex_pair> sedges;
  for (const auto * face : faces())
  {
    const auto fedges = face->edges();
    for (const auto & edge : fedges)
        sedges.insert(std::minmax(edge.first, edge.second));
  }
  return std::vector<vertex_pair> (sedges.begin(), sedges.end());
}

}  // end namespace
