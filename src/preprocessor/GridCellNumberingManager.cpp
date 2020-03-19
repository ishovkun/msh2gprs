#include "GridCellNumberingManager.hpp"

namespace gprs_data {

GridCellNumberingManager::GridCellNumberingManager(const mesh::Mesh & grid)
    : _grid(grid)
{}

DoFNumbering * GridCellNumberingManager::get_cell_numbering()
{
  size_t icell = 0;
  DoFNumbering * pdn = new DoFNumbering;
  pdn->m_cells.resize( _grid.n_cells(), pdn->m_unmarked );
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
    pdn->m_cells[cell->index()] = icell++;

  return pdn;
}


}  // end namespace gprs_data
