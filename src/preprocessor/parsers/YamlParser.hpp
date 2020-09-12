#pragma once
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep

#include "PreprocessorConfig.hpp"

#include <string>

namespace Parsers
{

class YamlParser
{
 public:
  YamlParser();
  void parse_file(const std::string & fname);
  PreprocessorConfig & get_config() {return config;}

 private:
  // sections
  void section_mesh(const YAML::Node & node);
  // var type can be 0 or 1: flow domain or geomechanics domain
  void section_domain_props(const YAML::Node &         node,
                            const ExpressionDomainType var_type);
  void embedded_fracs(const YAML::Node & node);
  void discrete_fracs(const YAML::Node & node);
  void boundary_conditions(const YAML::Node & node);
  void section_wells(const YAML::Node & node);
  void section_multiscale(const YAML::Node & node);
  // subsections
  void boundary_conditions_faces(const YAML::Node & node);
  void boundary_conditions_nodes(const YAML::Node & node);
  void subsection_cartesian_grid(const YAML::Node & node);
  // subsections
  void embedded_fracture(const YAML::Node       & node,
                         EmbeddedFractureConfig & conf);
  void discrete_fracture(const YAML::Node       & node,
                         DiscreteFractureConfig & conf);
  void domain(const YAML::Node &         node,
              const ExpressionDomainType var_type,
              DomainConfig     &         conf);
  void read_well(const YAML::Node & node,
                 WellConfig & conf);
  // subsubsection
  void bc_face(const YAML::Node & node, BCConfig & conf);
  void bc_node(const YAML::Node & node, BCConfig & conf);

  DomainConfig & get_domain_config(const int label);

  template <typename ValueType>
  void extract_subnode_value(const std::string           & key,
                             const  YAML::const_iterator & it,
                             ValueType                   & value)
  {
    try {
      value = it->second[key].as<ValueType>();
    }
    catch (YAML::TypedBadConversion<int> & error)
    {
      std::cout << "cannot convert " << key << ": "
                <<it->second[key] << std::endl;
      std::cout << error.what() << std::endl;
      std::cout << "aborting" << std::endl;
      exit(-1);
    }
  }

  // ATTRIBUTES
  PreprocessorConfig config;
};

}
