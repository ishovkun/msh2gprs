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
  // simple wells are vertical and only have one segment
  void setup_simple_well_(WellConfig const & conf);
  // complex wells are assigne with pairs of coordinates that indicate line segments
  // additionally, perforations might be specifically indicated
  void setup_complex_well_(WellConfig const & conf);
  // generate boundaries of the domain
  void generate_bounding_box_();
  double find_characteristic_length_() const;

  // config parameters for the grid generator
  INSIMMeshConfig const & _config;
  // well config
  std::vector<WellConfig> const & _wells;
  // this vector holds points to compute the bounding box for the grid
  std::vector<angem::Point<3,double>> _vertices;
  // holds mapping from well idx to its grid nodes
  std::unique_ptr<angem::Hexahedron<double>> _bbox{nullptr};
};

}  // end namespace gprs_data
