#pragma once
#include "config/WellConfig.hpp"
#include "mesh/Mesh.hpp"
#include "intersections/GridIntersectionSearcher.hpp"

namespace gprs_data {

/*
** This class searches for well nodes in the grid and saves them.
*/
class INSIMWellManager {
 public:
  INSIMWellManager(std::vector<WellConfig> const & wells, mesh::Mesh const & grid,
                   GridIntersectionSearcher const & searcher);

  // returns list of well vertices.
  // first component - well index. second component - grid vertex index.
  // single well can span over several grid vertices.
  std::vector<std::vector<size_t>> const & get_well_vertices() const {return _well_vertex_indices;}

  virtual ~INSIMWellManager() = default;

 private:
  void find_well_vertices_();

  std::vector<WellConfig> const & _wells;
  mesh::Mesh const & _grid;
  GridIntersectionSearcher const & _searcher;
  std::vector<std::vector<size_t>> _well_vertex_indices;
};

}  // end namespace gprs_data
