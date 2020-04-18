#pragma once

#include "mesh/Cell.hpp"
#include "mesh/Mesh.hpp"

namespace discretization {

/**
 * This class implements subdivisions/triangulatioins of mesh::Cell class.
 * Given a cell mesh::Cell entity, this class creates a new grid that is
 * built by dividing a cell into tetrahedral.
 * Keep in mind, that this class does not perform any grid refinement and
 * does not modify the original grid.
 * The method is exemplified (with pictures) in
 * Bishop, A displacement-based finite element formulation for
 * general polyhedra using harmonic shape functions (2014).
 *
 */
class Subdivision {
 public:
  /**
   * Constructor.
   * Takes a reference to a grid cells.
   * The order argument designates the number of recursive refinements
   * performed on the grid.
   * Input:
   * \param[in] cell : a grid cell to build triangulation for
   * \param[in] order : a number of recursive subdivisions for the triangulation
   * \param[out] triangulation : a grid formed by the triangulation and subsequent subdivisions
   *                             of cell.
   */
  Subdivision(const mesh::Cell & cell, mesh::Mesh & triangulation, const size_t order = 0);

 private:
  const mesh::Cell & _parent_cell;
  mesh::Mesh & _triangulation;
  const size_t _order;
};


}  // end namespace discretization
