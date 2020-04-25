#pragma once

namespace bc_functions {


/*  Returns true if |x-x_target| < tol */
double almost_equal(const double x, const double x_target, const double tol);

/*  Returns true if sqrt((x-x_target)^2 + (y - y_target)^2 + (z-z_target)^2) < tol */
double near(const double x, const double x_target,
            const double y, const double y_target,
            const double z, const double z_target, const double tol);

}  // end namespace gprs_data
