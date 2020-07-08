#pragma once
#include "FEMMethod.hpp"
#include "PolyhedralFEMSubdivision.hpp"
#include "SolverType.hpp"

enum PolyhedronIntegrationRule
{
  Full,
  FacesAverage,
  VerticesAverage,
  FacesPointwise,
  VerticesPointwise
};


struct FiniteElementConfig
{
  FEMMethod method = FEMMethod::strong_discontinuity;
  PolyhedralFEMSubdivision subdivision_method = PolyhedralFEMSubdivision::refinement;
  double solver_tolerance = 1e-5;                 // tolerance for msrsb convergence
  SolverType solver = SolverType::direct;
  size_t order = 0;  // refinement order (the more the finer)
  PolyhedronIntegrationRule integration_rule = PolyhedronIntegrationRule::Full;
};
