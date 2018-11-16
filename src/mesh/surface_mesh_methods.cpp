#include <surface_mesh_methods.hpp>

#include <cmath>      // std::pow

namespace mesh
{
const std::size_t MAX_EDGES = estimate_max_edges();

std::size_t estimate_max_edges()
{
  return std::sqrt(std::pow(2, 64) - 1);
}



}
