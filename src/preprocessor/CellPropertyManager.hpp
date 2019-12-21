#pragma once

#include "PreprocessorConfig.hpp"
#include "SimData.hpp"

namespace gprs_data {

class CellPropertyManager
{
 public:
  CellPropertyManager(const CellPropertyConfig & cell_properties,
                      SimData & data);

 private:
  const CellPropertyConfig & config;
  SimData & data;
};

}  // end namespace gprs_data
