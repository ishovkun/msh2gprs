#pragma once
#include "logger/Logger.hpp"
#include <memory>  // unique_ptr
#include <experimental/filesystem>

namespace gprs_data {

struct GlobalOpts
{
  using Path = std::experimental::filesystem::path;
  Path config_file_path;
  std::string log_file;
  logging::LogLevel log_level = logging::LogLevel::Debug;
  bool print_progressbar = true;

  static GlobalOpts & ref()
  {
    if (!_opts) _opts = std::unique_ptr<GlobalOpts>(new GlobalOpts());
    return *_opts;
  }

  // -------------- DO NOT TOUCH THESE --------------
 private:
  GlobalOpts() { _opts = nullptr; }
  static std::unique_ptr<GlobalOpts> _opts;
};

}  // end namespace gprs_data
