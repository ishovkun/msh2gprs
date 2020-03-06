#pragma once

#include <cstddef>  // provides std::size_t
#include <limits>  // provides std::numeric_limits

namespace mesh
{

namespace constants {

const std::size_t invalid_index = std::numeric_limits<std::size_t>::max();

// these marker are user for splitting cells with a plane
const int marker_below_splitting_plane = -1;
const int marker_above_splitting_plane = -2;
const int marker_splitting_plane = -3;
const int default_face_marker = -1;
const int default_cell_marker = -1;

}

}  // end namespace mesh
