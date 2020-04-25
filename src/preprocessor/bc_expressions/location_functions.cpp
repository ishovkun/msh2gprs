#include "location_functions.hpp"
#include <cmath>

namespace bc_functions {

double almost_equal(const double x, const double x_target, const double tol)
{
  return std::fabs(x - x_target) < tol;
}

/*  Returns true if sqrt((x-x_target)^2 + (y - y_target)^2 + (z-z_target)^2) < tol */
double near(const double x, const double x_target,
            const double y, const double y_target,
            const double z, const double z_target, const double tol)
{
  const double dx = x - x_target;
  const double dy = y - y_target;
  const double dz = z - z_target;
  return std::sqrt(dx*dx + dy*dy + dz*dz) < tol;
}

}  // end namespace gprs_data
