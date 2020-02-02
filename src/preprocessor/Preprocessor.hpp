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

class Preprocessor
{
 public:
  Preprocessor(const Path config_file_path);
  void run();

 private:
  void build_flow_discretization_();
  void combine_flow_discretizations_();
  void build_geomechanics_discretization_();
  void read_config_file_(const Path config_file_path);
  void read_mesh_file_(const Path mesh_file_path);
  void create_output_dir_();
  void write_output_();
  void build_dfem_discretization_();

  // ------------------ Variables ----------------------- //
  PreprocessorConfig config;
  SimData data;
  Path m_output_dir;
  std::shared_ptr<CellPropertyManager> pm_property_mgr;
  std::shared_ptr<DiscreteFractureManager> pm_dfm_mgr;
  std::shared_ptr<EmbeddedFractureManager> pm_edfm_mgr;
  std::shared_ptr<DoFNumbering> pm_flow_dof_numbering;
};

}  // end namespace gprs_data
