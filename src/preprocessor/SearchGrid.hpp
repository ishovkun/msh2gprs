#pragma once
#include "mesh/Mesh.hpp"
#include "angem/Hexahedron.hpp"

namespace gprs_data {

/*  Implements a search grid.
 *  Idea: We superimpose a cartesian uniform grid onto the original unstructured
 *  grid. We map each cell of the original grid to the superimosed grid.
 *  This structure is used for fast checking of collision detection.
 */
 class SearchGrid {
  public:
   SearchGrid(const mesh::Mesh & grid);
   size_t n_cells() const noexcept;
   std::array<int, 3> index(const angem::Point<3,double> & p) const;
   size_t global_index(const std::array<int, 3> &loc) const;
   size_t global_index(int i, int j, int k) const;
   std::array<int, 3> local_index(size_t idx) const;
   bool in_bounds(const angem::Point<3,double> & p) const noexcept;
   bool in_bounds(int i, int j, int k) const noexcept;
   std::vector<size_t> neighbors(size_t search_cell) const;

  private:
   void map_cell_(size_t cell_index);
   void compute_stepping_and_find_origin_();
   angem::Hexahedron<double> get_cell_hex_(size_t idx) const;
   void add_neighbor_(int i, int j, int k, std::vector<size_t> & dst) const;

   const mesh::Mesh & _grid;
   std::vector<std::list<size_t>> _mapping;
   angem::Point<3,double> _stepping;  // grid cell sizes
   angem::Point<3,double> _origin;
   std::array<int, 3> _dims;
 };


}  // end namespace gprs_data
