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

  FEPointData() {}

  FEPointData(size_t n_vertices)
      : values(n_vertices, 0.f), grads(n_vertices), weight(0.f) {}
};

struct FiniteElementData
{
  std::vector<FEPointData> points;  // gauss integration points
  FEPointData              center;  // values in cell center
  size_t                   element_index;
  std::shared_ptr<angem::Polyhedron<double>> topology;

  FiniteElementData() {}

  explicit FiniteElementData(size_t n_vertices, size_t n_points)
      : points(n_points, FEPointData(n_vertices)),
        center(n_vertices)
  {}
};

}  // end namespace discretization
