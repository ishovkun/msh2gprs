#include "FEMFaceDoFManager.hpp"

namespace discretization {

DoFNumbering FEMFaceDoFManager::build(const mesh::Mesh & grid, const std::vector<size_t> & face_indices) const
{
  DoFNumbering nmb;
  nmb.m_vertices.resize( grid.n_vertices(), nmb.m_unmarked );
  size_t dof = 0;
  for (const size_t iface : face_indices)
  {
    const mesh::Face & face = grid.face(iface);
    for (const size_t v : face.vertices())
      if ( nmb.m_vertices[v] == nmb.m_unmarked)
        nmb.m_vertices[v] = dof++;
  }
  nmb.m_n_dofs = dof;
  return nmb;
}

}  // end namespace discretization
