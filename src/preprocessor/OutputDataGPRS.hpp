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
  void save_flow_data_(const std::string cv_file, const std::string con_file,
                       const std::string trans_update_file) const;
  // save flow CV data
  void save_control_volume_data_(std::ofstream & out) const;
  // dave transmissibilities
  void save_trans_data_(std::ofstream & out) const;
  // dave transmissibility formulas for geomechanics update
  void save_trans_update_formulas_(std::ofstream & out) const;

  // save everything related to geomechancis discretization
  void save_geomechanics_data_() const;
  // save geometry data for geomechanics discretization
  void save_geometry_() const;
  // save ordered cell vertices and cell types for geomechanics
  void save_cell_geometry_(std::ofstream & out, const mesh::Mesh & grid) const;
  // save ordered face vertices and face types for geomechanics
  void save_face_geometry_(std::ofstream & out, const mesh::Mesh & grid) const;
  // save geomechanics keywords
  void save_geomechanics_keywords_() const;
  // save computed element data: (grad) shape functions, gauss weights, JxW
  void save_fem_data_() const;
  void save_embedded_fractures_(const std::string file_name) const;
  // save dirichlet and neumann boundary conditions for geomechanics
  void save_geomechanics_boundary_conditions_() const;
  void save_dirichlet_component_vertices(const size_t comp, const std::string comp_name,
                                         std::ofstream & out) const;
  void save_discrete_fracture_properties_(const std::string file_name) const;
  void saveWells(const std::string file_name) const;
  void saveFlowMultiScaleData(const std::string file_name);
  void saveMechMultiScaleData(const std::string file_name);

protected:
  const SimData & _data;
  const GPRSOutputConfig _config;
  mutable std::string _output_path;
};

}
