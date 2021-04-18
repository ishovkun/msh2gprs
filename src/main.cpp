#include "preprocessor/Preprocessor.hpp"
#include "preprocessor/GlobalOpts.hpp"
#include "preprocessor/parsers/CommandLineParser.hpp"
#include "logger/Logger.hpp"
#include <iostream>
#include <string>

int main(int argc, char *argv[])
{
  try {
    Parsers::CommandLineParser cmd(argc, argv);

    auto const & opts = gprs_data::GlobalOpts::ref();
    auto & logger = logging::Logger::ref();
    logger.set_file(opts.log_file);
    logger.set_verbosity(opts.log_level);

    // read stuff
    gprs_data::Preprocessor preprocessor(opts.config_file_path);
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
