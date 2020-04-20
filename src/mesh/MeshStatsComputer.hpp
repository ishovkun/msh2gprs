#pragma once
#include "Mesh.hpp"

namespace mesh {

/**
 * This class computes some simple mesh statistics.
 */
class MeshStatsComputer
{
 public:
  MeshStatsComputer(const mesh::Mesh & grid);

  double get_average_edge_length() const { return _average_edge_length; }

 protected:
  void compute_edge_charachteristics_();
  void compute_totals_();

  const Mesh & _grid;
  double _average_edge_length, _minimum_edge_length, _maximum_edge_length;
  double _total_edge_length;
  size_t _n_active_cells, _n_cells, _n_inactive_cells;
  size_t _n_active_faces, _n_faces, _n_inactive_faces;
};

}  // end namespace mesh
