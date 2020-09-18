#pragma once
#include "Mesh.hpp"

namespace mesh {

/* Performs grid refinement based on the edge aspect ratio */
class RefinementAspectRatio {
 public:
  RefinementAspectRatio(Mesh & mesh, const double aspect_ratio, const size_t max_level = 100);

 private:
  void find_problematic_cells_();
  bool check_cell_(const size_t cell_id) const;

  Mesh & _grid;
  const double _ratio;
  const size_t _max_level;
};


}  // end namespace mesh
