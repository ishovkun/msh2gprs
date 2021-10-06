#pragma once
#include "mesh/Mesh.hpp"                                  // provides mesh::Mesh
#include "config/MeshConfig.hpp"
#include "config/WellConfig.hpp"
#include "angem/Hexahedron.hpp"
#include <memory>  // unique_ptr

namespace gprs_data {

class GridGeneratorINSIM {
 public:
  GridGeneratorINSIM(INSIMMeshConfig const & config, std::vector<WellConfig> const & wells);
  operator mesh::Mesh() const;

  virtual ~GridGeneratorINSIM() = default;

 private:
  void setup_simple_well_(WellConfig const & conf);
  void generate_bounding_box_();
  double find_characteristic_length_() const;

  INSIMMeshConfig const & _config;
  std::vector<WellConfig> const & _wells;
  std::vector<angem::Point<3,double>> _vertices;
  std::unique_ptr<angem::Hexahedron<double>> _bbox{nullptr};
};

}  // end namespace gprs_data
