#include "DoFManager.hpp"

namespace gprs_data {

DoFManager::DoFManager(mesh::Mesh & grid,
             const std::vector<int> dfm_markers,
             const std::vector<int> edfm_markers)
    : m_grid(grid),
      m_set_dfm_markers(dfm_markers.begin(), dfm_markers.end()),
      m_set_edfm_markers(edfm_markers.begin(), edfm_markers.end())
{}

std::shared_ptr<DoFNumbering> DoFManager::distribute_dofs()
{
  std::shared_ptr<DoFNumbering> p_dofs = std::make_shared<DoFNumbering>();
  p_dofs->m_cells.resize(m_grid.n_cells_total());
  p_dofs->m_faces.reserve(m_grid.n_faces());

  size_t dof = 0;
  const auto & grid = m_grid;
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
    if (face->neighbors().size() == 2)  // some degenerate meshes remove face neighbors
      if (is_dfm_(face->marker()))
        p_dofs->m_faces[face->index()] = dof++;
  const size_t n_dfm = dof;
  std::cout << "Split dofs: " << n_dfm  << " dfm faces"<< std::endl;

  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
    if (face->neighbors().size() == 2)
      if (is_edfm_(face->marker()))
        p_dofs->m_faces[face->index()] = dof++;
  const size_t n_edfm = dof - n_dfm;
  std::cout << "Split dofs: " << n_edfm  << " edfm faces"<< std::endl;

  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
    p_dofs->m_cells[cell->index()] = dof++;
  const size_t n_cell_cv = dof - n_dfm - n_edfm;
  std::cout << "Split dofs: " << n_cell_cv  << " cell cvs"<< std::endl;

  // total number of degrees of freedom
  p_dofs->m_n_dofs = dof;
  return p_dofs;
}

std::shared_ptr<DoFNumbering> DoFManager::distribute_unsplit_dofs()
{
  std::shared_ptr<DoFNumbering> p_dofs = std::make_shared<DoFNumbering>();
  const auto & grid = m_grid;
  size_t dof = 0;

  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
    if (face->neighbors().size() == 2)  // some degenerate meshes remove face neighbors
      if (is_dfm_(face->marker()))
        p_dofs->m_faces[face->index()] = dof++;

  // edfm faces
  for (auto face = grid.begin_active_faces(); face != grid.end_active_faces(); ++face)
    if (face->neighbors().size() == 2)
      if (is_edfm_(face->marker()))
        p_dofs->m_faces[face->index()] = dof++;

  // reservoir cells
  p_dofs->m_cells.resize(grid.n_cells_total());
  // NOTE: raw iterator
  for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
    if (cell->parent() == *cell)  // skip refined cells here
    {
      p_dofs->m_cells[cell->index()] = dof;
      for (const auto p_cell : cell->ultimate_children())
        p_dofs->m_cells[p_cell->index()] = dof;
      dof++;
    }

  p_dofs->m_n_dofs = dof;
  return p_dofs;
}

void DoFManager::remap(std::vector<discretization::ControlVolumeData> & cv_data,
                       std::vector<discretization::ConnectionData>    & connection_data,
                       const DoFNumbering                             & old_dofs,
                       const DoFNumbering                             & new_dofs)
{
  // copy old data
  const auto cv_old = cv_data;
  cv_data.clear(); cv_data.resize( new_dofs.n_dofs() );
  // copy cvs
  for (const auto & cv : cv_old)
  {
    if (cv.type == discretization::ControlVolumeType::cell)
      cv_data[ new_dofs.cell_dof(cv.master) ] = cv;
    else // (cv.type == discretization::ControlVolumeType::face)
      cv_data[ new_dofs.face_dof(cv.master) ] = cv;
  }
  // modify connections
  for (auto &con : connection_data)
  {
    //  for (size_t & idof : con.elements)
    for (std::size_t i = 0; i < con.elements.size(); ++i)
    {
      if (cv_old[con.elements[i]].type == discretization::ControlVolumeType::cell)
        con.elements[i] = new_dofs.cell_dof(con.elements[i]);
      else // if( cv_old[con.elements[i]].type == discretization::ControlVolumeType::face )
        con.elements[i] = new_dofs.face_dof(con.elements[i]);
    }
  }
}

}  // end namespace gprs_data
