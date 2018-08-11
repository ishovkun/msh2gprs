#pragma once
#include <Polygon.hpp>

struct DomainConfig
{
  int model;
  int label;
  std::vector<std::string> expressions;
  std::vector<std::string> variables;
  // map expressions to variables
  std::map<int,int> local_to_global_vars;
  std::map<int,int> global_to_local_vars;

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
  // all variables used for function parsing
  std::vector<std::string> all_vars = {"x", "y", "z"};
  // output file names
  std::string domain_file = "domain.txt";
  std::string efrac_file = "efrac.txt";
  std::string bcond_file = "bcond.txt";
};


template<typename T>
std::size_t find(const T & item, const std::vector<T> & vec)
{
  for (std::size_t i=0; i<vec.size(); ++i)
    if (vec[i] == item)
      return i;
  return vec.size();
};
