#pragma once

#include <Shape.hpp>
#include <json.hpp>
#include <vector>

namespace Parsers
{


struct domain_config
{
  int model;
  double biot;
  double porosity;
  double permeability;
  double density;
  double young_modulus;
  double poisson_ratio;
  double thermal_expansion;
  double heat_capacity;
  double temperature;
  double pressure;
  double ref_temperature;
  double ref_pressure;
};


struct bc_config
{
  int type;
  angem::Point<3,double> value;
};


struct simdata_config
{
  std::vector<angem::Shape<double>> fractures;
  std::vector<domain_config> domains;
  std::vector<bc_config> bcs;
};


class Parser
{
 public:
  Parser();
  void parse_file(const std::string & fname);
  simdata_config & get_config();
 private:
  void parse_json(const std::string & fname);
  // ATTRIBUTES
  simdata_config config;
};

}
