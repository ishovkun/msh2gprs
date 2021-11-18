#pragma once
#include  "../config/MeshConfig.hpp"
#include "angem/Hexahedron.hpp"
#include <random>

namespace gprs_data {

/*
 * Implements Mitchell's best-candidate algorithm for adding "imaginary" wells.
 * The algorithm is outline in Appendix B. of
 * Guo, Zhenyu, and Albert C. Reynolds.
 * INSIM-FT in three-dimensions with gravity. Journal of Computational Physics 380 (2019): 143-169.
 */
class MitchellBestCandidate {
 public:
  // constructor
  MitchellBestCandidate(INSIMMeshConfig const & config);
  // main method
  [[nodiscard]] std::vector<angem::Point<3,double>> generate_points(std::vector<angem::Point<3,double>> const & actual,
                                                                    angem::Hexahedron<double> const & bounding_box);
  // destructor
  ~MitchellBestCandidate() = default;

 private:
  // uniformly sample points in 3d space so that lie within the bounding box
  // the number of points is taken as the size of ans
  // NOTE: I assume that the hexahedron aligns with xyz
  void sample_points_(angem::Hexahedron<double> const & bounding_box,
                      std::vector<angem::Point<3,double>> &ans);
  void compute_min_distances_(std::vector<angem::Point<3,double>> &candidates,
                              std::vector<angem::Point<3,double>> &actual,
                              std::vector<double> &distances,
                              size_t n_placed) const;

  INSIMMeshConfig const & _config;
  std::mt19937 _rnd;
};

}  // end namespace gprs_data
