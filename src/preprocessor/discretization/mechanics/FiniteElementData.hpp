#pragma once

#include "angem/Point.hpp"
#include "angem/Polyhedron.hpp"
#include <vector>
#include <memory>

namespace discretization {

struct FEPointData
{
  std::vector<double> values;
  std::vector<angem::Point<3,double>> grads;
  double weight;
};

struct FiniteElementData
{
  std::vector<FEPointData> points;  // gauss integration points
  FEPointData              center;  // values in cell center
  size_t                   element_index;
};

struct FiniteElementDataTopology
{
  FiniteElementData data;
  std::unique_ptr<angem::Polyhedron<double>> topology;
};

}  // end namespace discretization
