#pragma once
#include <Shape.hpp>

struct DomainConfig
{
  int model;
  int label;
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


struct SimdataConfig
{
  std::vector<angem::Shape<double>> fractures;
  std::vector<DomainConfig> domains;
  std::vector<BCConfig> bcs;
};
