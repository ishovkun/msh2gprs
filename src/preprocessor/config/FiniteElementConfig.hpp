#pragma once
#include "FEMMethod.hpp"
#include "PolyhedralFEMSubdivision.hpp"
#include "SolverType.hpp"


struct FiniteElementConfig
{
  FEMMethod method = FEMMethod::strong_discontinuity;
  PolyhedralFEMSubdivision subdivision_method = PolyhedralFEMSubdivision::refinement;
  double solver_tolerance = 1e-5;                 // tolerance for msrsb convergence
  SolverType solver = direct;
  size_t order;  // refinement order (the more the finer)
};
