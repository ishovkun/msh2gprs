#pragma once

#include "PreprocessorConfig.hpp"

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
  void parse_config_(const Path config_file_path);


  // ------------------ Variables ----------------------- //
  PreprocessorConfig m_config;
};

}  // end namespace gprs_data
