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
  void domain_props_json(const nlohmann::json::iterator & section_it);
  void embedded_fracs_json(const nlohmann::json::iterator & section_it);
  void boundary_conditions_json(const nlohmann::json::iterator & section_it);
  std::pair<std::string,std::string>
  get_pair_json(const nlohmann::json::iterator & section_it);
  // ATTRIBUTES
  SimdataConfig config;
  std::string comment = "_comment_";
};

}
