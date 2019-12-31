#include "active_cell_iterator.hpp"

namespace mesh
{

active_cell_iterator::active_cell_iterator(Cell * p_cell)
    : p_cell(p_cell)
{
  if (p_cell)
    assert (p_cell->m_parent != constants::invalid_index &&
            "Trying to get active cell iterator for inactive cell");
}

bool active_cell_iterator::operator==(const active_cell_iterator & other) const
{
  if (p_cell && other.p_cell)
    return *p_cell == *other.p_cell;
  else
    return p_cell == other.p_cell;
}

bool active_cell_iterator::operator!=(const active_cell_iterator & other) const
{
  return !operator==(other);
}

active_cell_iterator &
active_cell_iterator::operator++()
{
  while ( !p_cell->is_active() )
  {
    if (p_cell->index() == p_cell->m_grid_cells.size() - 1)
    {
      p_cell = nullptr;
      return *this;
    }
    const size_t next_cell_index = p_cell->index() + 1;
    p_cell = &(p_cell->m_grid_cells[next_cell_index]);
  }
  return *this;
}
 
}   // end namespace mesh
