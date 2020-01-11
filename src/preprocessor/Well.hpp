#pragma once

#include "angem/Point.hpp"
#include "PreprocessorConfig.hpp"

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
  bool simple() const {return segments.empty();}

  // user-input data
  std::vector<std::pair<angem::Point<3,double>, angem::Point<3,double>>> segments;
  // whether a segment is perforated
  std::vector<bool> perforated;
  angem::Point<3,double> coordinate;
  double radius, reference_depth;
  std::string name;
  bool reference_depth_set = false;

  // well data after computing mesh collision
  // flow volumes: cells and fractures intersected by the well
  std::vector<std::size_t> connected_volumes;
  // segment length at each intersection
  std::vector<double> segment_length;
  std::vector<angem::Point<3,double>> directions;
  // productivity at each intersection
  std::vector<double> indices;
};
