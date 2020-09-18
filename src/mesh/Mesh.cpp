#include <Mesh.hpp>
#include <SurfaceMesh.hpp>
#include "Cell.hpp"
#include "angem/PolyhedronFactory.hpp"
#include "angem/Collisions.hpp"    // angem::split
#include "yaml/include/yaml-cpp/emittermanip.h"
#include <unordered_set>
#include <algorithm>  // std::max
#include  <numeric>   // iota

namespace mesh
{

using std::vector;

Mesh::Mesh()
    : _n_inactive_cells(0)
{}

std::size_t Mesh::insert_cell(const std::vector<std::size_t> & ivertices,
                              const int                        vtk_id,
                              const int                        marker)
{
  // get faces: vectors of vertexindices based on vtk numbering
  const std::vector<std::vector<std::size_t>> poly_faces =
      angem::PolyhedronFactory::get_global_faces<double>(ivertices, vtk_id);
  std::vector<FaceTmpData> face_vector(poly_faces.size());
  for (size_t i = 0; i < poly_faces.size(); ++i)
  {
    face_vector[i].vertices = poly_faces[i];
    face_vector[i].vtk_id = face_vtk_id_(poly_faces[i].size());
  }

  std::vector<size_t> take_faces(face_vector.size());
  std::iota(take_faces.begin(), take_faces.end(), 0);
  return insert_cell_(ivertices, take_faces, face_vector,vtk_id, marker);
}

std::size_t Mesh::
insert_cell_(const std::vector<size_t> take_faces,
             const std::vector<FaceTmpData> &big_face_vector,
             const int                        marker)
{
  // cells don't have many vertices: don't need a hashmap here
  std::set<size_t> unique_vertices;
  for (const size_t iface : take_faces)
    for (const size_t v : big_face_vector[iface].vertices)
      unique_vertices.insert(v);

  std::vector<size_t> ivertices(unique_vertices.begin(), unique_vertices.end());
  return insert_cell_(ivertices, take_faces, big_face_vector,
                      angem::VTK_ID::GeneralPolyhedronID, marker);
}

std::size_t Mesh::
insert_cell_(const std::vector<std::size_t> & ivertices,
             const std::vector<size_t> take_faces,
             const std::vector<FaceTmpData> &big_face_vector,
             const int                        vtk_id,
             const int                        marker)
{
  // parent so that function can be reused for adding cells after splitting
  const std::size_t new_cell_index = n_cells_total();
  std::vector<std::size_t> face_indices;
  face_indices.reserve(take_faces.size());

  for (const size_t iface : take_faces)
  {
    const auto & face = big_face_vector[iface];
    const std::size_t face_index = insert_face_(face);
    face_indices.push_back(face_index);
  }

  m_cells.emplace_back(new_cell_index, ivertices, std::move(face_indices),
                       m_vertices, m_cells, m_faces, vtk_id, marker);

  for (const size_t vertex: ivertices)
  {
    m_vertex_cells[vertex].push_back(new_cell_index);
  }

  return new_cell_index;
}

size_t Mesh::insert_face(const std::vector<std::size_t> & ivertices,
                         const int                        vtk_id,
                         const int                        marker,
                         const std::size_t                face_parent)
{
  FaceTmpData f;
  f.vertices = ivertices;
  f.vtk_id = vtk_id;
  f.marker = marker;
  f.parent = face_parent;
  return insert_face_(f);
}

size_t Mesh::insert_face_(const FaceTmpData & f)
{
  if (m_vertex_cells.size() < n_vertices())
    m_vertex_cells.resize(n_vertices());
  if (m_vertex_faces.size() < n_vertices())
    m_vertex_faces.resize(n_vertices());

  size_t face_index = find_face(f.vertices);
  if (face_index == constants::invalid_index)
  {
    face_index = m_faces.size();
    m_faces.emplace_back(face_index, f.vertices, f.vtk_id, f.marker,
                         m_cells, m_faces, m_vertices, m_vertex_cells, f.parent);
    if (f.parent != constants::invalid_index)
    {
      assert( f.parent < n_faces() );
      m_faces[f.parent].m_children.push_back(face_index);
    }

  }
  else {
    if (f.parent != constants::invalid_index)
    {
      if (f.parent != face_index)  //  this sometimes happens when splitting cell through hanging nodes edge
      {
        m_faces[f.parent].m_children.push_back(face_index);
        m_faces[face_index].m_parent = f.parent;
      }
    }
    if (f.marker != constants::default_face_marker)
    {
      m_faces[face_index].m_marker = f.marker;
    }

    // assert(m_faces[face_index].m_vtk_id == f.vtk_id);
  }
  for (const size_t vertex : f.vertices)
  {
    validate_vertex_(vertex);
    if (!std::count( m_vertex_faces[vertex].begin(), m_vertex_faces[vertex].end(), face_index))
      m_vertex_faces[vertex].push_back(face_index);
  }

  return face_index;
}

active_cell_const_iterator Mesh::begin_active_cells() const
{
  for (auto cell = begin_cells(); cell != end_cells(); ++cell)
    if (cell->is_active()) return active_cell_const_iterator(&*cell);
  return active_cell_const_iterator(nullptr);
}

active_cell_iterator Mesh::begin_active_cells()
{
  for (auto cell = begin_cells(); cell != end_cells(); ++cell)
    if (cell->is_active()) return active_cell_iterator(&*cell);
  return active_cell_iterator(nullptr);
}

active_face_const_iterator Mesh::begin_active_faces() const
{
  for (auto face = begin_faces(); face != end_faces(); ++face)
    if (!face->is_active()) continue;
    else return active_face_const_iterator(&*face, m_faces);
  return active_face_const_iterator(nullptr, m_faces);  // end iterator
}

void Mesh::coarsen_cells()
{
  if ( n_cells_total() == n_active_cells() )  // no need to clear
    return;

  // find all deleted cells
  size_t min_cell_delete_index = std::numeric_limits<size_t>::max();
  for (auto cell = begin_cells(); cell != end_cells(); ++cell)
  {
    if (!cell->m_children.empty())
      cell->m_children.clear();
    else if (cell->parent() != *cell)  // to be deleted
      min_cell_delete_index = std::min(min_cell_delete_index, cell->index());
  }

  //  clear unused vertices
  for ( auto & vertex_cells : m_vertex_cells )
    for (auto it_cell = vertex_cells.begin(); it_cell != vertex_cells.end();)
      if (*it_cell >= min_cell_delete_index)
        vertex_cells.erase(it_cell);
      else ++it_cell;

  // find minimum vertex to erase: new vertices are always at the end
  size_t min_vertex_to_delete = std::numeric_limits<size_t>::max();
  for (std::size_t i=0; i<m_vertex_cells.size(); ++i)
    if (m_vertex_cells[i].empty())
    {
      min_vertex_to_delete = i;
      break;
    }
  m_vertex_cells.erase(m_vertex_cells.begin() + min_vertex_to_delete, m_vertex_cells.end() );

  // clear faces: if face has a deleted vertex then delete it
  // these faces are also consequtive and put into the end
  size_t min_face_to_delete = std::numeric_limits<size_t>::max();
  for (const auto & vertex_faces : m_vertex_faces)
    for (const size_t iface : vertex_faces )
      min_face_to_delete = std::min( min_face_to_delete, iface );
  m_faces.erase( m_faces.begin() + min_face_to_delete, m_faces.end() );
  m_vertex_faces.erase( m_vertex_faces.begin() + min_vertex_to_delete, m_vertex_faces.end() );
  // delete cells
  m_cells.erase( m_cells.begin() + min_cell_delete_index, m_cells.end() );
  // no more split cells
  _n_inactive_cells = 0;
}

std::pair<std::vector<size_t>, std::vector<size_t>>
separate_into_unique_groups(const std::vector<size_t> &group1,
                            const std::vector<size_t> &group2)
{
  std::pair<std::vector<size_t>, std::vector<size_t>> result;
  for (const size_t v : group1)
    if (!std::count(group2.begin(), group2.end(), v))
      result.first.push_back(v);
  for (const size_t v : group2)
    if (!std::count(group1.begin(), group1.end(), v))
      result.second.push_back(v);
  return result;
}

std::vector<size_t> find_vertices_from_both_groups(const std::vector<size_t> &verts,
                                                   const std::vector<size_t> &group1,
                                                   const std::vector<size_t> &group2)
{
  std::vector<size_t> result;
  for (const size_t v : verts)
  {
    if (std::count(group1.begin(), group1.end(), v))
    {
      result.push_back(v);
      break;
    }
  }

  for (const size_t v : verts)
  {
    if (std::count(group2.begin(), group2.end(), v))
    {
      result.push_back(v);
      break;
    }
  }
  return result;
}

size_t Mesh::find_face(const std::vector<size_t> & face_vertices) const
{
  const std::vector<size_t> sorted_vertices = sort_copy_(face_vertices);
  for (const size_t vertex: face_vertices)
  {
    // vertex hasn't been mapped to faces yet
    if (vertex >= m_vertex_faces.size())
      return constants::invalid_index;

    for (const size_t iface: m_vertex_faces[vertex])
    {
      assert( n_faces() > iface );
      const Face & f = m_faces[iface];
      const std::vector<size_t> face_vertices = sort_copy_(f.vertices());
      if ( std::equal( face_vertices.begin(), face_vertices.end(),
                       sorted_vertices.begin(), sorted_vertices.end()) )
      {
        return f.index();
      }
    }
  }
  return constants::invalid_index;
}

Mesh & Mesh::operator=(const Mesh & other)
{
  m_vertices = other.m_vertices;
  m_cells = other.m_cells;
  m_faces = other.m_faces;
  m_vertex_cells = other.m_vertex_cells;
  m_vertex_faces = other.m_vertex_faces;
  _n_inactive_cells = other._n_inactive_cells;
  for (auto & cell : m_cells)
  {
    cell.pm_grid_vertices = &m_vertices;
    cell.pm_grid_cells = &m_cells;
    cell.pm_grid_faces = &m_faces;
  }

  for (auto &face : m_faces)
  {
    face.pm_grid_cells = &m_cells;
    face.pm_grid_faces = &m_faces;
    face.pm_grid_vertices = &m_vertices;
    face.pm_grid_vertex_cells = & m_vertex_cells;
  }
  return *this;
}

std::vector<size_t> Mesh::neighbors_indices_(const vertex_pair & edge) const
{
  std::vector<size_t> ncs;
  for (const size_t v1_cell : m_vertex_cells[edge.first])
    if (cell(v1_cell).is_active())
      if (std::count(m_vertex_cells[edge.second].begin(), m_vertex_cells[edge.second].end(), v1_cell))
        ncs.push_back( v1_cell );
  return ncs;
}

size_t Mesh::insert_vertex(const angem::Point<3,double> & coord)
{
  m_vertices.push_back(coord);
  m_vertex_cells.resize(n_vertices());
  m_vertex_faces.resize(n_vertices());
  return n_vertices() - 1;
}

void Mesh::reserve_vertices(size_t nv)
{
  if (nv > n_vertices())
  {
    m_vertices.reserve(nv);
    m_vertex_cells.reserve(nv);
    m_vertex_faces.reserve(nv);
  }
}

void Mesh::validate_vertex_(size_t v) const
{
  if (v >= n_vertices())
    throw std::invalid_argument("invalid vertex index " + std::to_string(v));
}

std::vector<const Face*>  Mesh::vertex_faces(size_t vertex_index) const
{
  validate_vertex_(vertex_index);
  std::vector<const Face*> result;
  for (const size_t face_index : m_vertex_faces[vertex_index])
  {
    const auto & face = m_faces[face_index];
    if (face.is_active())
      result.push_back(&face);
  }
  return result;
}

}  // end namespace mesh
