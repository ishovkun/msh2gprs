#pragma once
#include "angem/Point.hpp"
#include <vector>
#include <string>

struct WellConfig
{
  std::string name;
  double radius;
  std::vector<angem::Point<3,double>> coordinates;
  std::vector<bool> perforated;
  // use this flag only to force connecting the well to fractures in 2D
  bool force_connect_fractures = false;
};
