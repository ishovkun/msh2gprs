#include "DoFManager.hpp"

namespace gprs_data {

DoFManager::DoFManager(SimData & data,
                       const std::vector<int> dfm_markers,
                       const std::vector<int> edfm_markers)
    : m_data(data),
      m_set_dfm_markers(dfm_markers.begin(), dfm_markers.end()),
      m_set_edfm_markers(edfm_markers.begin(), edfm_markers.end())
{}

DoFNumbering DoFManager::distribute_split_dofs()
{
  discretization::DoFNumbering dofs;
  dofs.m_cells.resize(m_data.grid.n_cells());
  dofs.m_faces.reserve(m_data.grid.n_faces());

  size_t dof = 0;
  const auto & grid = m_data.grid;
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
    if (is_dfm_(face->marker()))
      dofs.m_faces[face->index()] = dof++;
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
    if (is_edfm_(face->marker()))
      dofs.m_faces[face->index()] = dof++;
  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
    dofs.m_cells[cell->index()] = dof++;

  // total number of degrees of freedom
  dofs.m_n_dofs = dof;
  return dofs;
}

DoFNumbering DoFManager::distribute_unsplit_dofs()
{
  DoFNumbering dofs;
  const auto & grid = m_data.grid;
  size_t dof = 0;

  // dfm faces
  std::unordered_map<size_t, std::vector<size_t>> face_children;
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
    if (is_dfm_(face->marker()))
    {
     const size_t parent_index = grid.ultimate_parent(*face).index();
     face_children[parent_index].push_back( face->index() );
    }
  for (const auto & pair_parent_children : face_children)
  {
    for (const size_t child : pair_parent_children.second)
      dofs.m_faces[child] = dof;
    dof++;
  }

  // edfm faces
  for (const int marker : m_set_edfm_markers)
  {
    std::unordered_map<size_t,std::vector<size_t>> cell_frac_faces;
    for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
      if (face->marker() == marker)
      {
        const auto & neighbor_cell = *face->neighbors()[0];
        const auto & parent_cell = neighbor_cell.ultimate_parent();
        cell_frac_faces[parent_cell.index()].push_back(face->index());
      }

    for (const auto & frac_faces : cell_frac_faces)
    {
      for (const auto &iface : frac_faces.second)
        dofs.m_faces[iface] = dof;
      dof++;
    }
  }

  // reservoir cells
  dofs.m_cells.resize(grid.n_cells());
  // NOTE: raw iterator
  for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
    if (cell->parent() == *cell)  // skip refined cells here
    {
      for (const size_t icell : cell->ultimate_children())
        dofs.m_cells[icell] = dof;
      dof++;
    }

  dofs.m_n_dofs = dof;
  return dofs;
}

}  // end namespace gprs_data
