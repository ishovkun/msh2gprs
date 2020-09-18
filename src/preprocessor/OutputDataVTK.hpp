#pragma once

#include "SimData.hpp"
#include "mesh/Mesh.hpp"         // mesh::Mesh
#include "mesh/io/VTKWriter.hpp"  // VTKWriter

namespace gprs_data
{

/** 
 * Implements output into VTK format to visualize the exported data in Paraview. 
 * The vtk output data is also used by the postprocessor to visualize the output 
 * of AD-GPRS.
 **/
class OutputDataVTK
{
 public:
  /**
   * Constructor.
   * Saves references to SimData container and filenames to be output.
   * 
   * @param  {SimData} sim_data       : data container
   * @param  {VTKOutputConfig} config : container for relevant file names
   */
  OutputDataVTK(const SimData & sim_data, const VTKOutputConfig config);
  void write_output(const std::string & output_path) const;

 private:
  void save_reservoir_flow_data_(const std::string & fname) const;
  void save_reservoir_mechanics_data_(const std::string & fname) const;
  void save_dfm_data(const std::string & fname) const;
  void save_edfm_data(const std::string & fname) const;
  void save_wells_(const std::string & fname) const;
  // size is different for mech and flow
  // void saveMultiScaleSupport(const multiscale::MultiScaleOutputData & ms,
  //                            const std::size_t                        size,
  //                            const std::string                      & prefix,
  //                            std::ofstream                          & out);

  const SimData & _data;
  const mesh::Mesh & m_flow_grid;
  const mesh::Mesh & m_mech_grid;
  VTKOutputConfig m_config;
};

}
