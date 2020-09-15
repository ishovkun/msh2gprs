#pragma once

#include "PreprocessorConfig.hpp"
#include "SimData.hpp"
#include "discretization/flow/DoFNumbering.hpp"

namespace gprs_data {

/** This class manages discretization of wells
 * Given a config with well coordinates and geometrical
 * properties, it finds the cells intersected by wells
 * and computes J indices.
 */
class WellManager
{
 public:
  /**
   * Constructor 
   * 
   * @param  config : config for wells
   * @param  data   : container for output data and reservoir properties
   * @param  dofs   : flow degrees of freedom
   * @param  flow_method : method for modeling embedded fractures
   */
  WellManager(const std::vector<WellConfig> & config,
              SimData & data,
              const discretization::DoFNumbering & dofs,
              const EDFMMethod edfm_method);
  /**
   * Build wells discretization
   */
  void setup();

 protected:
  void setup_simple_well_(Well & well);
  void setup_simple_well_fast_(Well & well);
  void setup_segmented_well_(Well & well);
  void compute_well_index_(Well & well);
  std::array<double,3> get_bounding_box_(const std::size_t icell) const;
  bool setup_simple_well_matrix_(Well & well, size_t cell_index);
  // returns false if no intersection found
  void setup_simple_well_to_fracture_(Well & well, size_t cell_index);
  void compute_WI_matrix_(Well & well, discretization::WellSegment & segment);
  void compute_WI_frac_(Well & well, discretization::WellSegment & segment);

  const std::vector<WellConfig> _config;
  SimData & _data;
  const discretization::DoFNumbering & _dofs;
  const EDFMMethod _edfm_method;
  // std::vector<std::vector<size_t>> _well_connected_cells;
};

}  // end namespace preprocessor
