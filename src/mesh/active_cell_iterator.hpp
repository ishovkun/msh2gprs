#pragma once

#include "Cell.hpp"
#include "active_cell_const_iterator.hpp"

namespace mesh
{

/* this class implements iterations over the active
 * cells of the Mesh class.
 * A cell may become inactive if it has been split with the
 * Mesh::split_cell method. */
class active_cell_iterator
{
 public:
  // constructor
  active_cell_iterator(Cell * cell);
  // construct active cell iterator from const iterator
  active_cell_iterator(const active_cell_const_iterator& other);
  // comparison
  bool operator==(const active_cell_iterator & other) const;
  // comparison
  bool operator!=(const active_cell_iterator & other) const;
  // incrementing
  active_cell_iterator & operator++();
  // derementing
  // active_cell_iterator & operator--();
  // access operator
  inline Cell & operator*() { return *p_cell; }
  // access operator
  inline Cell * operator->() { return p_cell; }

 protected:
  void increment_raw_iterator_();
  Cell * p_cell;  // pointer to the cell
};

}   // end namespace mesh
