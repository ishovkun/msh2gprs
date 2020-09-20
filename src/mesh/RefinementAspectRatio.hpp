#pragma once
#include "Mesh.hpp"
#include "CellSplitter.hpp"

namespace mesh {

/* Performs grid refinement based on the edge aspect ratio */
class RefinementAspectRatio {
 public:
  RefinementAspectRatio(Mesh & mesh, const double aspect_ratio, const size_t max_level = 100);

 private:
  void split_cells_();
  bool check_cell_(const size_t cell_id) const;
  angem::Plane<double> find_cut_plane_(const size_t cell_id) const;
  angem::Point<3,double> edge_center_(const vertex_pair & edge) const;
  double edge_length_(const vertex_pair & edge) const;


  Mesh & _grid;
  const double _ratio;
  const size_t _max_level;
  CellSplitter _splitter;
};


}  // end namespace mesh
