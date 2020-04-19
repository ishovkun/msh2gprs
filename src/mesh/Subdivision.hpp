#pragma once

#include "Cell.hpp"
#include "Mesh.hpp"

namespace mesh {

/**
 * This class implements subdivisions/triangulatioins of mesh::Cell class.
 * Given a cell mesh::Cell entity, this class creates a new grid that is
 * built by dividing a cell into tetrahedral.
 * Keep in mind, that this class does not perform any grid refinement and
 * does not modify the original grid.
 * The method is exemplified (with pictures) in
 * Bishop, A displacement-based finite element formulation for
 * general polyhedra using harmonic shape functions (2014).
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
   * \param[out] triangulation : a grid formed by the triangulation and subsequent
   *                             subdivisions of cell.
   * \param[in] order : a number of recursive subdivisions for the triangulation
   */
  Subdivision(const Cell & cell, Mesh & triangulation, const size_t order = 0);

 private:
  void create_master_cell_();
  void perform_subdivision_r0_(Cell & cell);
  void perform_subdivision_tetra_(Cell & cell);
  void insert_tetra_(const std::vector<size_t> & local_vertex_indices,
                     const std::vector<size_t> & global_vertex_indices,
                     const std::vector<std::vector<size_t>> & faces,
                     const std::vector<size_t> & face_parents,
                     const size_t parent_cell_index);
  std::vector<size_t> get_face_order_(const Cell & cell,
                                      const std::vector<size_t> &cell_vertices) const;

  // *********************** Veriables ****************************** //
  const Cell & _parent_cell;
  Mesh & _grid;
  const size_t _order;
  angem::PointSet<3, double> _created_vertices;
  std::vector<size_t> _created_vertex_indices;
};


}  // end namespace discretization
