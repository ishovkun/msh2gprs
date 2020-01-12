#pragma once

#include "PreprocessorConfig.hpp"
#include "SimData.hpp"

#include <experimental/filesystem>  // filesystem

namespace gprs_data {

namespace filesystem = std::experimental::filesystem;
using Path = filesystem::path;

class Preprocessor
{
 public:
  Preprocessor(const Path config_file_path);
  void run();

 private:
  void read_config_file_(const Path config_file_path);
  void read_mesh_file_(const Path mesh_file_path);
  void create_output_dir_();
  void write_output_();

  // ------------------ Variables ----------------------- //
  PreprocessorConfig config;
  SimData data;
  Path m_output_dir;
};

}  // end namespace gprs_data
