#pragma once
#include "angem/Point.hpp"
#include <array>
#include <cstddef>  // size_t
#include <limits>   // numeric_limits

namespace discretization {

/*  This structure stores a single perforation segment:
 *  a connection between a well and a flow control volume.
 */
struct WellSegment
{
  // static constexpr size_t dof_undefined = std::numeric_limits<size_t>::max();
  size_t dof;  // flow CV number the segment is connected to
  size_t element_id;  // id of connected cell or face
  std::array<double,3> bounding_box;  // temporary structure used to store bounding box dimensions
  angem::Point<3,double> direction;   // direction of the wekk segment
  double length = 0;              // length of the perforation
  bool perforated = true;
  double wi = 0;  // well index
};

}  // end namespace discretization
