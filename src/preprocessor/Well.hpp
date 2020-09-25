#pragma once

#include "angem/Point.hpp"
#include "PreprocessorConfig.hpp"
#include "discretization/wells/WellSegment.hpp"

#include <utility>  // pair
#include <vector>
#include <string>

/* Class that controls well geometry and productivity indices
 *
 * There are two types of wells: simple and complex
 * Simple wells are vertical and only have (x, y) coordinates
 * Complex wells are defined with pairs of (x, y, z) segments
  */
class Well
{
 public:
  Well(const WellConfig & config);
  bool simple() const noexcept {return segments.empty();}
  bool force_fracture_connection() const noexcept {return _force_frac_connect;}

  // user-input data
  std::vector<std::pair<angem::Point<3,double>, angem::Point<3,double>>> segments;
  // whether a segment is perforated
  std::vector<bool> perforated;
  angem::Point<3,double> coordinate;
  double radius, reference_depth;
  std::string name;
  bool reference_depth_set = false;

  bool _force_frac_connect;
  std::vector<discretization::WellSegment> segment_data;
};
