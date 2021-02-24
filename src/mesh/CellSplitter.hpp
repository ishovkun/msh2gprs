#pragma once

#include "Mesh.hpp"

namespace mesh {

/**
 * This class implements splitting of cells for cEDFM method.
 * It splits specified cells with planes by dividing them in two.
 */
class CellSplitter {
 public:

  /* Constructor:
   * \param[in] grid : a reference to the grid object to be split
   * \param[in] tol  : relative geometric tolerance to distinguish new split vertices
   */
  CellSplitter(Mesh & grid, double tol = 1e-4);

  /* Split a cell by cutting it with a plane. New cell indices are appended
   * at the back, so it is safe to split multiple cells in a row.
   *
   * Returns true if the split occurred
   *
   * Note: cell is copied since inserting new cells invalidates the pointers. */
  bool split_cell(Cell cell, const angem::Plane<double> & plane,
                  const int splitting_face_marker = constants::marker_splitting_plane);

  /* Set sink for debug output */
  void set_debug_output(std::ostream * os) {_os = os;}

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

  void track_new_vertices_(std::vector<size_t> & global_vertex_indices,
                           angem::PolyGroup<double> & split,
                           std::vector<size_t> & new_vertices);

  void create_face_groups_(angem::PolyGroup<double> & split,
                           const std::vector<std::vector<size_t>> & face_vertex_global_numbering,
                           const std::vector<size_t> & polygroup_polygon_parents,
                           const std::vector<Face*> & cell_faces,
                           const int splitting_face_marker,
                           const size_t split_face_local_index,
                           std::vector<size_t> & cell_above_faces,
                           std::vector<size_t> & cell_below_faces,
                           std::vector<FaceTmpData> & tmp_faces);

  Mesh & _grid;
  double const _tol;
  angem::PointSet<3, double> _new_vertex_coord;
  std::vector<size_t> _new_vertices;
  std::ostream * _os{nullptr};  // pointer to debug output stream
};

}  // end namespace mesh
