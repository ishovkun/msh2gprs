#pragma once
#include "mesh/Mesh.hpp"
#include "SimData.hpp"

namespace multiscale {

class MSFlow {
 public:
  MSFlow(mesh::Mesh const & grid, gprs_data::SimData & data);
  virtual ~MSFlow() = default;

 protected:
  void find_coarse_center(size_t coarse, std::vector<size_t> & cells);
  mesh::Mesh const & _grid;
  gprs_data::SimData & _data;
};

}  // end namespace multiscale
