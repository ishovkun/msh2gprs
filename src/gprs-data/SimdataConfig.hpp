#pragma once

#include <Polygon.hpp>

#include <map>
#include <memory> // shared / unique_ptr

struct DomainConfig
{
  int label;
  std::vector<std::string> expressions;
  std::vector<std::string> variables;
  // map expressions to variables
  std::map<int,int> local_to_global_vars;
  std::map<int,int> global_to_local_vars;
  bool coupled = true;
};


struct BCNodeConfig
{
  static const int type = 1;
  angem::Point<3,double> value;
  angem::Point<3,double> coord;
};


struct BCConfig
{
  int label;
  int type;
  angem::Point<3,double> value;
};


struct EmbeddedFractureConfig
{
  std::shared_ptr<angem::Polygon<double>> body;  // embedded fractures
  std::size_t n1 = 0, n2 = 0;                    // remeshing params
  double cohesion = 0;
  double friction_angle = 30;
  double dilation_angle = 0;
  double aperture = 1e-3;
  double conductivity = 10;
};


struct DiscreteFractureConfig
{
  int label;
  double conductivity = 10;
  double aperture = 1e-3;
};


struct WellConfig
{
  std::string name;
  double radius;
  std::vector<angem::Point<3,double>> coordinates;
  std::vector<bool> perforated;
};


struct SimdataConfig
{
  std::vector<EmbeddedFractureConfig>  fractures;  //  embedded  fractures
  std::vector<DiscreteFractureConfig>  discrete_fractures;
  std::vector<DomainConfig>            domains;
  std::vector<BCConfig>                bc_faces;
  std::vector<BCNodeConfig>            bc_nodes;
  std::vector<WellConfig>              wells;

  // all variables used for function parsing
  std::vector<std::string> all_vars = {"x", "y", "z"};
  std::vector<int>         expression_type;  // (0 - flow, 1 - mechanics, -1 -no output)
  std::vector<std::string> special_keywords =
  {"PERM", "PERMX", "PERMY", "PERMZ", "PORO", "VFACTOR"};
  static constexpr double default_permeability = 1;
  static constexpr double default_volume_factor = 1;

  // output file names
  std::string mesh_file;
  std::string domain_file           = "domain.txt";
  std::string mechanics_domain_file = "gm_geometry.txt";
  std::string efrac_file            = "gm_SDA.txt";
  std::string discrete_frac_file    = "gm_DFM.txt";
  std::string bcond_file            = "bcond.txt";
  std::string wells_file            = "wells.txt";
  // special keywords needed for computing fluid data
  // (they are not outputted)
  static constexpr double nan = -999.999;
  double node_search_tolerance = 1e-10;
  double frac_cell_elinination_factor = 0.2;
};


// find an item in a vector
template<typename T>
std::size_t find(const T & item, const std::vector<T> & vec)
{
  for (std::size_t i=0; i<vec.size(); ++i)
    if (vec[i] == item)
      return i;
  return vec.size();
}


// more generic implement of the previous func
template<typename iterable>
std::size_t find(const typename iterable::value_type & item,
                 const iterable & container)
{
  std::size_t counter = 0;
  for (const auto & iter_item : container)
    if (iter_item == item)
      return counter;
    else
      counter++;
  return counter;
}
