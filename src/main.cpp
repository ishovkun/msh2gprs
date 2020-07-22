#include "preprocessor/Preprocessor.hpp"
#include <iostream>
#include <string>
#include "logger/Logger.hpp"
#include <experimental/filesystem>

namespace filesystem = std::experimental::filesystem;
using Path = filesystem::path;

int main(int argc, char *argv[])
{
  // process cmd arguments
  if (argc < 2)
  {
    std::cout << "please specify a config file."
              << std::endl
              << "Example: "
              << "msh2gprs config.json"
              << std::endl
              << "an example config is distributed with this code"
              << std::endl;
    return 0;
  }
  if (argc > 2)
  {
    std::cout << "Please provide only a single input argument" << std::endl;
    return 1;
  }

  // config file
  const std::string fname_config = argv[1];
  const Path path_config(fname_config);

  auto & logger = logging::Logger::ref();
  const Path log_file_path = path_config.parent_path() / "log.txt";
  logger.set_file(log_file_path.filename());
  logger.set_verbosity(logging::LogLevel::Debug);

  try {
    // read stuff
    gprs_data::Preprocessor preprocessor(path_config);
    // let the fun begin
    preprocessor.run();
  }
  catch (const std::exception& e)
  {
    logging::critical() << e.what() << std::endl;
    return 1;
  }

  return 0;
}
