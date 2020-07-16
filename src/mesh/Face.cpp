#include "Face.hpp"
#include "Cell.hpp"

namespace mesh
{

Face::Face(const std::size_t                       face_index,
           const std::vector<std::size_t> &        face_vertices,
           const int                               face_vtk_id,
           const int                               face_marker,
           std::vector<Cell>              &        grid_cells,
           std::vector<Face>                     & grid_faces,
           std::vector<Point> &                    grid_vertices,
           std::vector<std::vector<std::size_t>> & grid_vertex_cells,
           const std::size_t                       parent)
    : m_index(face_index),
      m_vertices(face_vertices),
      m_vtk_id(face_vtk_id),
      m_marker(face_marker),
      pm_grid_cells(&grid_cells),
      pm_grid_faces(&grid_faces),
      pm_grid_vertices(&grid_vertices),
      pm_grid_vertex_cells(&grid_vertex_cells),
      m_parent(parent)
{
  if (m_parent == constants::invalid_index)
    m_parent = index();
}


Face & Face::operator=(const Face & other)
{
  m_index = other.index();
  m_vertices = other.vertices();
  m_vtk_id = other.vtk_id();
  m_marker = other.marker();
  pm_grid_cells = other.pm_grid_cells;
  pm_grid_faces = other.pm_grid_faces;
  pm_grid_vertices = other.pm_grid_vertices;
  pm_grid_vertex_cells = other.pm_grid_vertex_cells;
  m_parent = other.parent().index();
  return *this;
}

std::vector<const Cell*> Face::neighbors() const
{
  std::vector<const Cell*> result;
  std::vector<const Cell*> const_neibs = raw_neighbors();
  for (const Cell* p_cell : const_neibs)
  {
    if (p_cell->is_active())
      result.push_back(p_cell);
  }

  if (result.size() > 2 || result.empty())
  {
    std::cout << std::endl;
    std::cout << "face index " << index() << std::endl;
    std::cout << "face vertices: ";
    for (auto v : vertices())
      std::cout << v << " ";
    std::cout << std::endl;
    for (const Cell* p_cell : const_neibs)
      std::cout << "tried " << p_cell->index() << std::endl;
    std::cout << "is_active() = " << is_active() << std::endl;
    std::cout << "neibs ";
    for (const auto p_cell : result)
      std::cout << p_cell->index() << " ";
    std::cout << " (";
    for (const auto p_cell : result)
      std::cout << p_cell->ultimate_parent().index() << " ";
    std::cout << " )";
    std::cout << std::endl << std::flush;
  }
  assert(!result.empty());
  assert(result.size() <= 2);
  // sort result by index
  // This is done to avoid changing the order of neighbors
  // if vertex reordering is performed
  std::sort(result.begin(), result.end(), [](const auto c1, const auto c2) -> bool
                                          {
                                            return c1->index() < c2->index();
                                          });
  return result;
}

std::vector<Cell*> Face::neighbors()
{
  std::vector<Cell*> result;
  const auto & cthis = *this;
  std::vector<const Cell*> const_neibs = cthis.neighbors();
  for (const Cell* p_cell : const_neibs)
    result.push_back( const_cast<Cell*>(p_cell) );

  return result;
}

std::vector<const Cell*> Face::raw_neighbors() const
{
  // find how many times each vertex neighboring cell is encountered
  std::unordered_map<size_t, size_t> cell_time;
  for (const size_t vertex : m_vertices)
  {
    for (const size_t cell_index: (*pm_grid_vertex_cells)[vertex])
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
    {
      const Cell & cell = (*pm_grid_cells)[it.first];
      face_neighbors.push_back(&cell);
    }
  }

  return face_neighbors;
}

std::vector<Point> Face::vertex_coordinates() const
{
  std::vector<Point> coordinates;
  coordinates.reserve(m_vertices.size());
  for (const std::size_t vertex_index : vertices())
    coordinates.push_back( (*pm_grid_vertices)[vertex_index] );
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

double Face::area() const
{
  return polygon().area();
}


Point Face::center() const
{
  return polygon().center();
}

bool Face::has_vertex(const std::size_t vertex_index) const
{
  if (std::find(m_vertices.begin(), m_vertices.end(), vertex_index) ==
      m_vertices.end())
    return false;
  else return true;
}

std::vector<vertex_pair> Face::edges() const
{
  const auto & verts = vertices();
  std::vector<vertex_pair> pairs(verts.size());
  for (std::size_t i=0; i<verts.size(); ++i)
  {
    std::size_t i1, i2;
    if (i < verts.size() - 1)
    {
      i1 = verts[i]; i2 = verts[i+1];
    }
    else
    {
      i1 = verts[ i ]; i2 = verts[ 0 ];
    }
    pairs[i] = std::make_pair(i1, i2);
  }
  return pairs;
}

const Face & Face::ultimate_parent() const
{
  const Face * par =  &parent();
  while (par->parent() != *par)
    par = &par->parent();
  return *par;
}

Face & Face::ultimate_parent()
{
  const auto & cthis = *this;
  return const_cast<Face&>(cthis.ultimate_parent());
}

bool Face::has_edge(const vertex_pair & edge) const
{
  const vertex_pair edge_sorted = std::minmax(edge.first, edge.second);
  for (const auto e : edges())
  {
    const vertex_pair es = std::minmax(e.first, e.second);
    if ( edge_sorted.first == es.first && edge_sorted.second == es.second )
      return true;
  }
  return false;
}

}
