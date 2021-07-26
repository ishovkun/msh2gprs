#pragma once
#include "SimData.hpp"
#include "mesh/Mesh.hpp"
#include "intersections/UniformCartesianGrid.hpp"
#include <array>
#include <cstddef>

namespace multiscale {

class GeometricPartition {
 public:
  GeometricPartition(std::array<size_t,3> dims, gprs_data::SimData & data);
  ~GeometricPartition() = default;
  std::vector<size_t> get() const {return _partition;}

 private:
  gprs_data::UniformCartesianGrid build_uniform_grid_(std::array<size_t,3> dims) const;

  gprs_data::SimData const & _data;
  gprs_data::UniformCartesianGrid _uniform;
  std::vector<size_t> _partition;
};

}  // end namespace multiscale
