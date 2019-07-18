#pragma once

#include "angem/Shape.hpp"
#include <json.hpp>
#include <SimdataConfig.hpp>
#include <vector>

namespace Parsers
{

class JsonParser
{
 public:
  JsonParser();
  void parse_file(const std::string & fname);
  SimdataConfig & get_config();
 private:
  void parse(const std::string & fname);
  // var type can be 0 or 1: flow domain or geomechanics domain
  void domain_props(const nlohmann::json::iterator & section_it,
                         const int                        var_type);
  void embedded_fracs(const nlohmann::json::iterator & section_it);
  void boundary_conditions(const nlohmann::json::iterator & section_it);
  void boundary_conditions_faces(nlohmann::json::iterator it,
                                 const nlohmann::json::iterator & end);
  void boundary_conditions_nodes(nlohmann::json::iterator it,
                                 const nlohmann::json::iterator & end);
  void boundary_conditions_face(nlohmann::json::iterator it,
                                const nlohmann::json::iterator & end,
                                BCConfig & conf);
  void boundary_conditions_node(nlohmann::json::iterator it,
                                const nlohmann::json::iterator & end,
                                BCNodeConfig & conf);
  void embedded_fracture(nlohmann::json::iterator it,
                         const nlohmann::json::iterator & end,
                         EmbeddedFractureConfig & conf);
  void discrete_fracs(const nlohmann::json::iterator & section_it);
  void discrete_fracture(nlohmann::json::iterator it,
                         const nlohmann::json::iterator & end,
                         DiscreteFractureConfig & conf);



  std::pair<std::string,std::string>
  get_pair(const nlohmann::json::iterator & section_it);
  // ATTRIBUTES
  SimdataConfig config;
  std::string comment = "_comment_";
};

}
