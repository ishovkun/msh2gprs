#pragma once
#include "config/WellConfig.hpp"
#include "mesh/Mesh.hpp"
#include "intersections/GridIntersectionSearcher.hpp"
#include "discretization/flow/DoFNumbering.hpp"
#include "discretization/flow/ControlVolumeData.hpp"
#include "Well.hpp"

namespace gprs_data {

/*
** This class searches for well nodes in the grid and saves them.
*/
class INSIMWellManager {
 public:
  /*
  ** Constructor.
  ** Given configuration of wells, grid, and a grid searcher, build the list of vertices
  ** that lie closest to the well segment centers
   */
  INSIMWellManager(std::vector<WellConfig> const & wells, mesh::Mesh const & grid,
                   GridIntersectionSearcher & searcher);

  // returns list of wells.
  std::vector<Well> const & get_wells() const;
  // setup wells for ouput
  void assign_dofs(discretization::DoFNumbering const & dofs);
  //
  void compute_well_indices(std::vector<discretization::ControlVolumeData> const & cvs);
  // returns list of well vertices.
  // first component - well index. second component - grid vertex index.
  // single well can span over several grid vertices.
  std::vector<std::vector<size_t>> get_well_vertices() const;

  virtual ~INSIMWellManager() = default;

 private:
  std::vector<size_t> find_well_vertices_(Well const & well);
  void create_well_perforations_(Well & well, std::vector<size_t> const & well_vertices);
  void compute_bounding_box_(discretization::WellSegment & segment);


  void setup_wells_();

  std::vector<WellConfig> const & _config;
  mesh::Mesh const & _grid;
  GridIntersectionSearcher & _searcher;
  std::vector<std::vector<size_t>> _well_vertex_indices;
  std::vector<Well> _wells;
};

}  // end namespace gprs_data
