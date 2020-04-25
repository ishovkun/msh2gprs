#pragma once

namespace bc_functions {

/**
 * Consider a prismatic beam of length L and rectangular cross-section of width 2a and height 2b,
 * with coordinate system shown in Figure 11(a). The beam is fixed (weakly) on the end x3 D L and
 * subjected to a transverse shear force F in the negative x2-direction at the opposite end x3 D 0.
 * This function returns the analytical solution for x-displacement ux.
 */
double cantilever_beam_end_shear_ux(const double x, const double y, const double z,
                                    const double F, const double E, const double nu,
                                    const double a, const double b);
 /* This function returns the analytical solution for x-displacement uy. */
double cantilever_beam_end_shear_uy(const double x, const double y, const double z,
                                    const double F, const double E, const double nu,
                                    const double a, const double b);
 /* This function returns the analytical solution for x-displacement uz. */
double cantilever_beam_end_shear_uz(const double x, const double y, const double z,
                                    const double F, const double E, const double nu,
                                    const double a, const double b);

}  // end namespace gprs_data
