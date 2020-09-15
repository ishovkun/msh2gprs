#include "GridEntityNumberingManager.hpp"

namespace gprs_data {

GridEntityNumberingManager::GridEntityNumberingManager(const mesh::Mesh & grid)
    : _grid(grid)
{}

DoFNumbering * GridEntityNumberingManager::get_numbering()
{
  DoFNumbering * pdn = new DoFNumbering;
  // cell numbering
  pdn->m_cells.resize( _grid.n_cells_total(), pdn->m_unmarked );
  size_t icell = 0;
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell, ++icell)
    pdn->m_cells[cell->index()] = icell;

  // face numbering
  // pdn->m_faces.resize( _grid.n_faces(), pdn->m_unmarked );
  pdn->m_faces.reserve( _grid.n_faces() );
  size_t iface = 0;
  for (auto face = _grid.begin_active_faces(); face != _grid.end_active_faces(); ++face, ++iface)
    pdn->m_faces[face->index()] = iface;

  return pdn;
}


}  // end namespace gprs_data
