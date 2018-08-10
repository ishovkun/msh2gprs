#pragma once

#include <Shape.hpp>
#include <json.hpp>
#include <vector>
#include <SimdataConfig.hpp>

namespace Parsers
{

class Parser
{
 public:
  Parser();
  void parse_file(const std::string & fname);
  SimdataConfig & get_config();
 private:
  // JSON
  void parse_json(const std::string & fname);
  void rock_props_json(const nlohmann::json::iterator & section_it);
  void embedded_fracs_json(const nlohmann::json::iterator & section_it);
  void boundary_conditions_json(const nlohmann::json::iterator & section_it);
  // ATTRIBUTES
  SimdataConfig config;
};

}
