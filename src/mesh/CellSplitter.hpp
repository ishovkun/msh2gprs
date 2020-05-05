#pragma once

#include "Mesh.hpp"

namespace mesh {

/**
 * This class implements splitting of cells for cEDFM method.
 * It splits specified cells with planes by dividing them in two.
 */
class CellSplitter {
 public:
  /* Constructor */
  CellSplitter(Mesh & grid);
  /* Split a cell by cutting it with a plane. New cell indices are appended
   * at the back, so it is safe to split multiple cells in a row.
   * Note: cell is copied since inserting new cells invalidates the pointers. */
  void split_cell(Cell cell, const angem::Plane<double> & plane,
                  const int splitting_face_marker = constants::marker_splitting_plane);

 protected:
  /* get a vector of polygon global vertex indices given a vector with
   * local polygon vertex indices and a mapping vector. */
  std::vector<std::size_t> build_global_face_indices_(const std::vector<size_t> & polygon_local_indices,
                                                      const std::vector<size_t> & local_to_global) const;
  std::map<vertex_pair,size_t> find_affected_edges_(const std::vector<size_t> &new_vertices,
                                                    const Cell & cell) const;
  void insert_hanging_node_(const Cell parent, const vertex_pair edge,
                            const size_t inserted_vertex);
  /* if any of the cell faces contain the vertices in edge, split that face in two */
  void split_face_in_cell_(const Cell parent, const vertex_pair new_edge);


  Mesh & _grid;
  angem::PointSet<3, double> _new_vertex_coord;
  std::vector<size_t> _new_vertices;
};

}  // end namespace mesh
