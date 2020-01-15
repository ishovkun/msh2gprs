#pragma once

#include "PreprocessorConfig.hpp"
#include "SimData.hpp"
// #include "transes.hpp"

namespace gprs_data
{

class OutputDataGPRS
{
public:
  OutputDataGPRS(const SimData & data, const GPRSOutputConfig config);
  void write_output(const std::string & output_path) const;

 private:
  void save_flow_data_(const std::string cv_file, const std::string con_file) const;
  void saveGeometry(const std::string & output_path);
  void saveGeomechDataNewKeywords(const std::string file_name);
  void saveEmbeddedFractureProperties(const std::string file_name);
  void saveBoundaryConditions(const std::string file_name);
  void saveDiscreteFractureProperties(const std::string file_name);
  void saveWells(const std::string file_name) const;
  void saveFlowMultiScaleData(const std::string file_name);
  void saveMechMultiScaleData(const std::string file_name);

protected:
  const SimData & m_data;
  const GPRSOutputConfig m_config;
};

}
