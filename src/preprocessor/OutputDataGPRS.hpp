#pragma once

#include "PreprocessorConfig.hpp" // provide GPRSOutputConfig
#include "SimData.hpp"            // provide SimData

namespace gprs_data
{

/** This class implements an output writer for AD-GPRS format.
 * It takes data from SimData container and writes it into a bunch of files
 * in the appropriate format.
 **/
class OutputDataGPRS
{
public:
  /**
   * Constructor.
   * 
   * @param[in]  data   : container for data to be saved
   * @param[in]  config : file names for GPRS output
   */
  OutputDataGPRS(const SimData & data, const GPRSOutputConfig config);
  /**
   * Save data into output_path director.
   * @param  {std::string} output_path : string that describes the output path
   */
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
