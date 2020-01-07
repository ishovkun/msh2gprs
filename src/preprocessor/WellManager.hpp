#pragma once

#include "PreprocessorConfig.hpp"
#include "SimData.hpp"

namespace gprs_data {

class WellManager
{
 public:
  WellManager(const std::vector<WellConfig> & config,
              SimData & data);
  void setup();

 protected:
  void setup_simple_well_(Well & well);
  void setup_segmented_well_(Well & well);
  void compute_well_index_(Well & well);
  angem::Point<3,double> get_dx_dy_dz_(const std::size_t icell) const;

 private:
  const std::vector<WellConfig> m_config;
  SimData & m_data;
  std::vector<std::vector<size_t>> m_well_connected_cells;
};

}  // end namespace preprocessor
