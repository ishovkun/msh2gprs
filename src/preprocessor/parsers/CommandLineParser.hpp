#pragma once
#include "cxxopts.hpp"
#include "logger/Logger.hpp"

namespace Parsers
{

class CommandLineParser {
 public:
  CommandLineParser(int argc, char *argv[]);
  virtual ~CommandLineParser() = default;

 private:
  void parse_paths_();
  void parse_logger_();

  cxxopts::Options _options;
  cxxopts::ParseResult _result;
};

}
