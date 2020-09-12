#pragma once

#include "PreprocessorConfig.hpp"
#include "SimData.hpp"
#include "CellPropertyManager.hpp"
#include "DiscreteFractureManager.hpp"
#include "EmbeddedFractureManager.hpp"
#include "discretization/flow/DoFNumbering.hpp"

#include <experimental/filesystem>  // filesystem

namespace gprs_data {

namespace filesystem = std::experimental::filesystem;
using Path = filesystem::path;
using discretization::DoFNumbering;

/* This is the main class of the program.
 * It call various methods to read config file, grid, and then
 * do the preprocessing:
 * - split the grid with dfm fractures,
 * - distribute properties in the domain
 * - build flow and geomechanics discretization
 * - distribute wells and compute J indices
 * - compute multiscale data. */
class Preprocessor
{
public:
  /* Constructor.
   * During the construction, the config file and mesh file are read.
   * Takes a path to the configuration file as an argument. */
  Preprocessor(const Path config_file_path);
  /* Run the preprocessor. All the action is taking place here */ 
  void run();

private:
  // build flow discretization
  void build_flow_discretization_();
  // build geomechanics discretization
  void build_geomechanics_discretization_();
  // read yaml config file
  void read_config_file_(const Path config_file_path);
  // create grid
  void setup_grid_(const Path config_dir_path);
  // read .msh or .vtk grid file
  void read_mesh_file_(const Path mesh_file_path);
  // create an output directory. If it already exists, cleans it first.
  // WARNING: this can destroy your files.
  void create_output_dir_();
  // write output of the program
  void write_output_();

  // ------------------ Variables ----------------------- //
  PreprocessorConfig config; // config of the preprocessor written in yaml file
  SimData data; // data container. stores data for output.
  Path m_output_dir; // path to output directory
  std::shared_ptr<CellPropertyManager> pm_property_mgr;  // pointer to cell property manager
  std::shared_ptr<DiscreteFractureManager> pm_dfm_mgr;   // pointer to discrete fracture manager
  std::shared_ptr<DiscreteFractureManager> pm_cedfm_mgr;   // pointer to discrete fracture manager
  std::shared_ptr<EmbeddedFractureManager> pm_edfm_mgr;  // pointer to embedded fracture manager
};

}  // end namespace gprs_data
