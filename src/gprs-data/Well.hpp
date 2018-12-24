#pragma once

#include <angem/Point.hpp>
#include <SimdataConfig.hpp>

#include <utility>  // pair
#include <vector>
#include <string>

class Well
{
 public:
  Well(const WellConfig & config);
  bool simple() const {return segments.empty();}

  // attributes
  std::vector<std::pair<angem::Point<3,double>, angem::Point<3,double>>> segments;
  std::vector<bool> perforated;
  angem::Point<3,double> coordinate;
  double radius;
  std::string name;
  std::vector<std::size_t> connected_volumes;
};
