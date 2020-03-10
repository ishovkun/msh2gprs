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
   */
  WellManager(const std::vector<WellConfig> & config,
              SimData & data,
              const discretization::DoFNumbering & dofs);
  /**
   * Build wells discretization
   */
  void setup();

 protected:
  void setup_simple_well_(Well & well);
  void setup_segmented_well_(Well & well);
  void compute_well_index_(Well & well);
  angem::Point<3,double> get_bounding_box_(const std::size_t icell) const;

 private:
  const std::vector<WellConfig> m_config;
  SimData & m_data;
  const discretization::DoFNumbering & m_dofs;
  std::vector<std::vector<size_t>> m_well_connected_cells;
};

}  // end namespace preprocessor
