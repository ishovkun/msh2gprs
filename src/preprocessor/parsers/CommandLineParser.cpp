#include "CommandLineParser.hpp"
#include "../GlobalOpts.hpp"
#include <string>
#include <vector>
#include <iostream>

namespace Parsers
{

using namespace cxxopts;
using Path = std::experimental::filesystem::path;
using std::string;
using std::vector;

CommandLineParser::CommandLineParser(int argc, char *argv[])
    : _options("ERSP", "Easy Reservoir Simulation Preprocessor")
{
  _options
      .set_tab_expansion()
      .positional_help("/path/to/input/file.yaml")
      .add_options()
      // available input arguments
      ("i,input", "Input config file path", value<string>())
      //
      ("d,debug", "Log level debug", value<bool>()->default_value("true")) // a bool parameter
      //
      ("e,experimental", "Use features that are currently under development",
       value<bool>()->default_value("false"))
      //
      ("o,optimized", "Use optimized yet unreliable algorithms",
       value<bool>()->default_value("false"))
      //
      ("no_progressbar", "do not draw progressbars", value<bool>()->default_value("false"))
  ;

  _options.parse_positional({"input"});

  _result =_options.parse(argc, argv);

  try
  {
    parse_paths_();
    parse_logger_();
    gprs_data::GlobalOpts::ref().enable_experimental = _result["e"].as<bool>();
  }
  catch (option_has_no_value_exception const & e)
  {
    std::cout << _options.help() << std::endl;
    throw e;
  }
  catch (std::runtime_error const & e)
  {
    std::cout << _options.help() << std::endl;
    throw e;
  }

}

void CommandLineParser::parse_paths_()
{
  auto & opts = gprs_data::GlobalOpts::ref();
  if (!_result.count("input"))
    throw std::runtime_error("User must provide input file");

  const std::string fname_config = _result["i"].as<string>();
  opts.config_file_path = std::experimental::filesystem::path(fname_config);
  opts.log_file = (opts.config_file_path.parent_path() / "log.txt").filename();
}

void CommandLineParser::parse_logger_()
{
  auto & opts = gprs_data::GlobalOpts::ref();
  bool const debug_log = _result["d"].as<bool>();
  if (debug_log) opts.log_level = logging::LogLevel::Debug;
  opts.print_progressbar = !_result["no_progressbar"].as<bool>();
}


}
