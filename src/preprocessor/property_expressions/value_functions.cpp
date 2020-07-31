#include "value_functions.hpp"
#include <algorithm>  // min, max

namespace property_functions {

double clip(const double x, const double xmin, const double xmax)
{
  return std::min(std::max(x, xmin), xmax);
}

}  // end namespace property_functions
