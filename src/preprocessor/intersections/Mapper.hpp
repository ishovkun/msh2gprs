#pragma nce
#include "mesh/Mesh.hpp"
#include "UniformCartesianGrid.hpp"

namespace gprs_data {

/*  Implements a search grid.
 *  Idea: We superimpose a cartesian uniform grid onto the original unstructured
 *  grid. We map each cell of the original grid to the superimosed grid.
 *  This structure is used for fast checking of collision detection.
 */
 class Mapper {
  public:
   Mapper(const mesh::Mesh & grid);
   std::array<int, 3> index(const angem::Point<3,double> & p) const;
   void map_cell(size_t cell_index);

  private:
   UniformCartesianGrid build_grid_() const;

   const mesh::Mesh & _grid;
   UniformCartesianGrid _cartesian;
   std::vector<std::list<size_t>> _mapping;
 };


}  // end namespace gprs_data
