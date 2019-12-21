#include "Preprocessor.hpp"

#include <string>
#include "parsers/YamlParser.hpp"

namespace gprs_data {

Preprocessor::Preprocessor(const Path config_file_path)
{
  if (!filesystem::exists(config_file_path))
  {
    const std::string error_msg =
        "config file does not exist: " +
        std::string(filesystem::absolute(config_file_path));
    throw std::invalid_argument(error_msg);
  }

  parse_config_(config_file_path);
}

void Preprocessor::parse_config_(const Path config_file_path)
{
  const std::string fname = config_file_path.filename();
  const std::size_t str_len = fname.size();
  if (fname.substr(str_len - 4, str_len) == "yaml")
  {
    Parsers::YamlParser parser;
    
  }
}

}  // end namespace gprs_data
