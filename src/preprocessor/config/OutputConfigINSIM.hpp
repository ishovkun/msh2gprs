#pragma once
#include <string>

struct OutputConfigINSIM
{
  std::string cv_file    = "cell_data.txt";
  std::string tube_file  = "tube_data.txt";
  std::string wells_file = "wells.txt";
  double transmissibility_mult = 1.f;  // insim model takes transmissibility in md
};
