#pragma once
#include "angem/Point.hpp"
#include <array>
#include <cstddef>  // size_t

namespace discretization {

/*  This structure stores a single perforation segment:
 *  a connection between a well and a flow control volume.
 */
struct WellSegment
{
  size_t element;  // flow CV number the segment is connected to
  std::array<double,3> bounding_box;  // temporary structure used to store bounding box dimensions
  angem::Point<3,double> direction;   // direction of the wekk segment
  double segment_length;              // length of the perforation
  bool perforated;
};

}  // end namespace discretization
