#pragma once

#include <cassert>       // provides assert
#include <cstddef>       //std::size_t
#include <vector>        //std::vector
#include <unordered_map> //unordered_map

namespace gprs_data { class DoFManager; }

namespace discretization {

class DoFNumbering
{
 public:
  DoFNumbering() : m_n_dofs(0) {};
  inline size_t cell_dof(const size_t cell_index) const {
    assert( cell_index < m_cells.size() );
    return m_cells[cell_index];
  }

  inline size_t face_dof(const size_t face_index) const
  {
    const auto it_face = m_faces.find(face_index);
    assert (it_face != m_faces.end());
    return it_face->second;
  }

  inline size_t vertex_dof(const size_t vertex_index) const {
    assert( vertex_index < m_vertices.size() && "invalid vertex index" );
    return m_vertices[vertex_index];
  }

  inline bool has_vertex(const size_t vertex_index) const {
    if ( vertex_index >= m_vertices.size() )
      return false;
    else if ( m_vertices[vertex_index] ==  m_unmarked)
      return false;
    else return true;
  }

  inline bool is_active_face(const size_t face_index) const {return m_faces.find(face_index)!=m_faces.end();}
  inline size_t n_dofs() const { return m_n_dofs; }

 protected:
  std::vector<size_t> m_vertices;            // vertex dofs
  std::vector<size_t> m_cells;               // cell to cv
  std::unordered_map<size_t,size_t> m_faces; // face to cv
  size_t m_n_dofs;                           // total number of degrees of freedom
  const size_t m_unmarked = std::numeric_limits<size_t>::max();  // label for unmarked vertices

  // buds
  friend class gprs_data::DoFManager;
  friend class FEMFaceDoFManager;
};

}  // end namespace discretization
