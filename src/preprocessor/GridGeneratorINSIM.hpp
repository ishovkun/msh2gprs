#pragma once
#include "mesh/Mesh.hpp"                                  // provides mesh::Mesh

namespace gprs_data {

class GridGeneratorINSIM {
 public:
  GridGeneratorINSIM();
  operator mesh::Mesh() const;

  virtual ~GridGeneratorINSIM() = default;
};

}  // end namespace gprs_data
