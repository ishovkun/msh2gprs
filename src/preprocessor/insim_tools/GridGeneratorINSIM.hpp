#pragma once
#include "mesh/Mesh.hpp"                                  // provides mesh::Mesh
#include "config/MeshConfig.hpp"
#include "config/WellConfig.hpp"
#include "angem/Hexahedron.hpp"
#include <memory>  // unique_ptr

namespace gprs_data {

/* Implemnets grid generator more or less as discussed by Guo [2019].
 * Builds a bounding box for the point cloud of well vertices.
 * Adds padding to this bounding box.
 * Takes the center of each well segment and inserts it as am embedded vertex in the grid.
 * Generates a grid (currently by calling GMsh).
 */
class GridGeneratorINSIM {
 public:
  /* Constructor.
   * Input:
   * \param[in] config : configuration parameters for the grid generator
   * \param[in] wells  : geometric data for multiple wells
   */
  GridGeneratorINSIM(INSIMMeshConfig const & config, std::vector<WellConfig> const & wells);
  // generates and returns mesh generated from well data
  operator mesh::Mesh() const;

  virtual ~GridGeneratorINSIM() = default;

 private:
  // simple wells are vertical and only have one segment
  void setup_simple_well_(WellConfig const & conf);
  // complex wells are assigne with pairs of coordinates that indicate line segments
  // additionally, perforations might be specifically indicated
  void setup_complex_well_(WellConfig const & conf);
  // generate boundaries of the domain
  angem::Hexahedron<double> generate_bounding_box_() const;
  // add padding to bounding box so that wells are away from the domain boundries
  angem::Hexahedron<double> pad_bounding_box_(angem::Hexahedron<double> const &) const;
  // make the bounding box uniform so that GMsh takes a large grid size
  angem::Hexahedron<double> extend_bounding_box_(angem::Hexahedron<double> const &) const;
  // assign cell markers to the value specified in the config
  void assign_cell_labels_(mesh::Mesh & grid) const;
  // invoke Mitchell's algorithm
  // generate point candidates within the bounding box only
  void add_imaginary_wells_(angem::Hexahedron<double> const & box);

  // config parameters for the grid generator
  INSIMMeshConfig const & _config;
  // well config
  std::vector<WellConfig> const & _wells;
  // this vector holds points to compute the bounding box for the grid
  std::vector<angem::Point<3,double>> _vertices;
  // holds mapping from well idx to its grid nodes
  std::unique_ptr<angem::Hexahedron<double>> _bounding_box{nullptr};
};

}  // end namespace gprs_data
