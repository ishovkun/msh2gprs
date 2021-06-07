#pragma once
#include "mesh/Mesh.hpp"
#include "SimData.hpp"

namespace multiscale {

class MSFlow {
 public:
  MSFlow(mesh::Mesh const & grid, gprs_data::SimData & data, MultiscaleConfig const &config);
  virtual ~MSFlow() = default;

 protected:
  void find_coarse_center(size_t coarse, std::vector<size_t> & cells);
  std::function<double(double)> build_weight_function() const;
  mesh::Mesh const & _grid;
  gprs_data::SimData & _data;
  MSPartitioning const _type;
  size_t const _ncoarse;
};

}  // end namespace multiscale
