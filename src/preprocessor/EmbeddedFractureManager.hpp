#pragma once

#include "PreprocessorConfig.hpp" // provides EmbeddedFractureConfig, EDFMMethod
#include "discretization/flow/DoFNumbering.hpp" // provides DoFNumbering
#include "mesh/CellSplitter.hpp"
#include "SimData.hpp"  // provides SimData

namespace gprs_data {

/**
 * This class is responsible for Embedded fractures.
 * I can split the grid cells to create a cEDFM approach,
 * create a DFM-type configuration for embedded fracture to feed to discretization classes,
 * create a grid for vtk output,
 * map flow control volumes to edfm elements,
 * and map geomechanics cells to edfm elements.
 */
class EmbeddedFractureManager
{
 public:
  /**
   * Constructor.
   * A reference to SimData is made to be filled upon calling various functions.
   * 
   * @param[in]  config           : configuration for the class
   * @param[in]  edfm_method      : which edfm method to use
   * @param[in]  min_dist_to_node : minimum distance between fracture and grid vertex relative to cell size
   * @param[out] data             : container for output data
   */
  EmbeddedFractureManager(std::vector<EmbeddedFractureConfig> &config,
                          const EDFMMethod edfm_method,
                          const double min_dist_to_node,
                          SimData & data);

  /**
   * Split cells in the grid by embedded fractures.
   * This first computes which cells are affected and then asks
   * mesh::Mesh to split those cells.
   * The splits are performed in a sequence:
   * first, all cells for embedded fracture #1 are split,
   * second, all cells for embedded fracture #2 are split and so on.
   * Be careful with this method since it invalidates all properties on the grid.
   */
  void split_cells();
  /**
   * This function produces a vector of configurations for each embeddedd fracture
   * as if they were discrete fractures.
   * This method must be called after split_cells().
   * Use this function to feed it to discretization and obtain a 
   * cEDFM-type discretization.
   * 
   * @return {std::vector<DiscreteFractureConfig>}  : vector of configs for split embedded fractures
   */
  std::vector<DiscreteFractureConfig> generate_dfm_config();
  /**
   * Check whether a face with a marker face_marker is an embedded fracture
   * @param  {int} face_marker : face marker
   * @return {bool}            : true if face marker belongs to an edfm fracture after splitting cells
   */
  bool is_fracture(const int face_marker) const;
  /**
   * Generate geomechanical properties for embedded fracture (SDA formulation).
   * This fills out the sda_data parameter in SimData.
   */
  void distribute_mechanical_properties();
  /**
   * Get split face embedded fracture markers
   * NOTE: call this function only after split_cells()
   * @return {std::vector<int>}  : vector of unique face markers that describe split edfm cells
   */
  std::vector<int> get_face_markers() const;
  /**
   * Build edfm surface grid for vtk output
   * This fills SimData::edfm_grid and SimData::edfm_cell_mapping
   * @param  {discretization::DoFNumbering} dofs : flow degrees of freedom of edfm cells
   */
  // void build_edfm_grid(const discretization::DoFNumbering & dofs);
  // map SDA cells to edfm control volumes
  // do it only after coarsening the grid
  // and distribute mechanical properties
  void map_mechanics_to_control_volumes(const discretization::DoFNumbering & dofs);

 private:
  bool find_edfm_cells_(angem::Polygon<double> & fracture, std::vector<size_t> & cells);
  bool find_edfm_cells_fast_(angem::Polygon<double> & fracture, std::vector<size_t> & cells);
  void find_edfm_cells_and_faces_();
  // split internal grid cells due to intersection with embedded fracture
  void split_cells_(angem::Polygon<double> & fracture, std::vector<size_t> & cells, const int face_marker);
  // find the maximum face marker of the grid
  int find_maximum_face_marker_() const;
  // wrapper around m_marker_config;
  size_t fracture_index_(const int face_marker) const;
  // ------------------ Variables ----------------- //
  // non-const cause we move fractures to avoid collision with vertices
  std::vector<EmbeddedFractureConfig> &config;
  EDFMMethod m_method;                  // regular edfm, pedfm, or cedfm
  // minimum distance from fracture to vertex relative to the cell size
  const double _min_dist_to_node;
  SimData & m_data;                     // container for output data
  mesh::Mesh & m_grid;                  // reference to grid object
  std::map<int,size_t> m_marker_config; // marker to config index
  mesh::CellSplitter _splitter;
};

}  // end namespace gprs_data
