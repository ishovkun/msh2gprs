#include "active_cell_iterator.hpp"

namespace mesh
{

active_cell_iterator::active_cell_iterator(Cell * p_cell)
    : p_cell(p_cell)
{
  if (p_cell)
    assert (p_cell->is_active() &&
            "Trying to get active cell iterator for inactive cell");
}

active_cell_iterator::active_cell_iterator(const active_cell_const_iterator& other)
    : p_cell(&const_cast<Cell&>(*other))
{}

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
  if (p_cell->index() == p_cell->pm_grid_cells->size() - 1)
  {
    p_cell = nullptr;
    return *this;
  }

  increment_raw_iterator_();
  while ( !p_cell->is_active() )
  {
    if (p_cell->index() == p_cell->pm_grid_cells->size() - 1)
    {
      p_cell = nullptr;
      return *this;
    }
    increment_raw_iterator_();
  }
  return *this;
}

void active_cell_iterator::increment_raw_iterator_()
{
  const size_t next_cell_index = p_cell->index() + 1;
  p_cell = &((*(p_cell->pm_grid_cells))[next_cell_index]);
}

}   // end namespace mesh
