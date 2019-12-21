#pragma once

#include "PreprocessorConfig.hpp"
#include "SimData.hpp"

namespace gprs_data {

class EmbeddedFractureManager
{
 public:
  EmbeddedFractureManager(const std::vector<EmbeddedFractureConfig> &config,
                          const EDFMMethod edfm_method,
                          SimData & data);
  void split_cells();

 private:
  bool find_edfm_cells_(angem::Polygon<double> & fracture,
                        std::vector<size_t> & cells) const;
  // split internal grid cells due to intersection with
  // embedded fracture
  void split_cells_(angem::Polygon<double> & fracture,
                    std::vector<size_t> & cells);
  // ------------------ Variables -----------------
  const std::vector<EmbeddedFractureConfig> &config;
  // simple or pedfm
  EDFMMethod m_method;
  // will be filled
  SimData & data;
  mesh::Mesh m_split_grid;
};

}  // end namespace gprs_data
