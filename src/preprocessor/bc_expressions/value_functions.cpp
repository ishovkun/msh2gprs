#include "value_functions.hpp"
#include <limits>  // std::numeric_limits
#include <cmath>
#ifndef M_PI
#define M_PI  3.14159265358979323846  /* π */
#endif /* M_PI */

namespace bc_functions {

double cantilever_beam_end_shear_ux(const double x, const double y, const double z,
                                    const double F, const double E, const double nu,
                                    const double a, const double b)
{
  const double I = 4 * a * b*b*b / 3;
  return -F * nu / (E*I) * x * y * z;
}

double cantilever_beam_end_shear_uy(const double x, const double y, const double z,
                                    const double F, const double E, const double nu,
                                    const double a, const double b)
{
  const double I = 4 * a * b*b*b / 3;
  return F / (E*I) * (nu/2 * (x*x - y*y) * z - (z*z*z) / 6);
}

double minus_one_power(const std::size_t n)
{
  return (n % 2 == 0) ? 1 : -1;
}

double cantilever_beam_end_shear_uz(const double x, const double y, const double z,
                                    const double F, const double E, const double nu,
                                    const double a, const double b)
{
  const double I = 4 * a * b*b*b / 3;
  const double term1 = 0.5 * y * (nu * x*x + z*z);
  const double term2 = nu * (y*y*y) / 6;
  const double term3 = (1+nu) * (b*b * y - (y*y*y)/3);
  const double term4 = -(a*a) * nu * y / 3;
  const double term5m = -4 * (a*a*a) * nu / (M_PI*M_PI*M_PI);
  double summ = 0;
  for (std::size_t n=1; n<10; ++n)
    summ += minus_one_power(n) / (n*n*n) * cos(n * M_PI * x /a) *
                                 sinh(n*M_PI * y / a) / cosh(n*M_PI * y / a);
  return F/(E*I) * (term1 + term2 + term3 + term4 + term5m*summ);
}

}  // end namespace gprs_data
