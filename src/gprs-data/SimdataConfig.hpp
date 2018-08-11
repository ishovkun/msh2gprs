#pragma once
#include <Polygon.hpp>
#include <muParser.h>

struct DomainConfig
{
  int model;
  int label;
  mu::Parser function_parser;
  std::vector<std::string> expressions;
  // std::map<std::string, >
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


struct BCConfig
{
  int label;
  int type;
  angem::Point<3,double> value;
  const double nan = -999.999;
};

struct EmbeddedFractureConfig
{
  angem::Polygon<double> body;  // embedded fractures
  double cohesion = 0;
  double friction_angle = 30;
  double dilation_angle = 0;
};


struct SimdataConfig
{
  std::vector<EmbeddedFractureConfig> fractures;  // embedded fractures
  std::vector<DomainConfig> domains;
  std::vector<BCConfig> bcs;
};
