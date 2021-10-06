#pragma once
#include "mesh/Mesh.hpp"                                  // provides mesh::Mesh
#include "config/WellConfig.hpp"

namespace gprs_data {

class GridGeneratorINSIM {
 public:
  GridGeneratorINSIM(std::vector<WellConfig> const & wells);
  operator mesh::Mesh() const;

  virtual ~GridGeneratorINSIM() = default;

 private:
  std::vector<WellConfig> const & _wells;
};

}  // end namespace gprs_data
