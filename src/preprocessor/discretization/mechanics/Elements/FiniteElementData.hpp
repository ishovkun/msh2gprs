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

  FEPointData(size_t n_vertices) {resize(n_vertices);}

  void resize(size_t n_vertices)
  {
    values.assign(n_vertices, 0.f);
    grads.assign(n_vertices, {0.f, 0.f, 0.f});
    weight = 0.f;
  }
};

struct FiniteElementData
{
  std::vector<FEPointData> points;  // gauss integration points
  FEPointData              center;  // values in cell center
  size_t                   element_index;

  FiniteElementData() {}

  explicit FiniteElementData(size_t n_vertices, size_t n_points)
  { resize(n_vertices, n_points); }

  void resize(size_t n_vertices, size_t n_points)
  {
    points.assign(n_points, FEPointData(n_vertices));
    center.resize(n_vertices);
  }
};

}  // end namespace discretization
