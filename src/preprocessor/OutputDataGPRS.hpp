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
  // save data related to flow in the reservoir (but not wells)
  void save_flow_data_(const std::string cv_file, const std::string con_file) const;
  // save everything related to geomechancis discretization
  void save_geomechanics_data_() const;
  // save geometry data for geomechanics discretization
  void save_geometry_() const;
  // save ordered cell vertices for geomechanics
  void save_cell_geometry_(std::ofstream & out, const mesh::Mesh & grid) const;
  // save geomechanics keywords
  void save_geomechanics_keywords_() const;
  void saveEmbeddedFractureProperties(const std::string file_name);
  void saveBoundaryConditions(const std::string file_name);
  void saveDiscreteFractureProperties(const std::string file_name);
  void saveWells(const std::string file_name) const;
  void saveFlowMultiScaleData(const std::string file_name);
  void saveMechMultiScaleData(const std::string file_name);

protected:
  const SimData & _data;
  const GPRSOutputConfig _config;
  mutable std::string _output_path;
};

}
