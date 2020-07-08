#pragma once

enum class FEMMethod
{
  strong_discontinuity,      // just computes regular FEM shape functions
  polyhedral_finite_element, // Bishop, A displacement-based finite element formulation
                             // for general polyhedra using harmonic shape functions (2014).
  mixed                      // computes classic FEM shape funtions for regular cells and
                             // polyhedral shape functions for non-standard polyhedra
};
