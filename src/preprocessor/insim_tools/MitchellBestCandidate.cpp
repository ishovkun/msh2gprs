#include "MitchellBestCandidate.hpp"
#include <limits>

namespace gprs_data {

using Point = angem::Point<3,double>;

MitchellBestCandidate::MitchellBestCandidate(INSIMMeshConfig const & config)
    : _config(config),
      _rnd( std::random_device{}() )
{
  // enforce a consistent seed: good for debugging + we don't really need true randomness
  _rnd.seed(1);
}

std::vector<angem::Point<3,double>>
MitchellBestCandidate::generate_points(std::vector<angem::Point<3,double>> const & actual,
                                       angem::Hexahedron<double> const & box)
{
  std::vector<Point> ans( actual.size() + _config.n_imaginary_wells );
  std::copy( actual.begin(), actual.end(), ans.begin() );
  std::vector<Point> candidates( _config.n_candidates );
  std::vector<double> distances( _config.n_candidates );
  for (size_t i = 0; i < _config.n_imaginary_wells; ++i) {
    sample_points_( box, candidates );
    compute_min_distances_( candidates, ans, distances, actual.size() + i );
    auto const best = std::max_element( distances.begin(), distances.end() );
    size_t const best_idx = std::distance( distances.begin(), best );
    assert( best_idx < distances.size() );
    assert( actual.size() + i < ans.size() );
    ans[actual.size() + i] = candidates[best_idx];
  }
  return ans;
}

void MitchellBestCandidate::sample_points_(angem::Hexahedron<double> const & bounding_box,
                                           std::vector<angem::Point<3,double>> &ans)
{
  auto const & vertices = bounding_box.get_points();
  Point o = vertices[0];
  std::array<Point, 3> directions {
    vertices[1] - o,
    vertices[3] - o,
    vertices[4] - o,
  };

  std::uniform_real_distribution<double> distribution(0.f, 1.f);
  for (size_t i = 0; i < ans.size(); ++i) {
    Point p = o;
    for (size_t j = 0; j < 3; ++j)
      p += distribution(_rnd) * directions[j];
    ans[i] = p;
  }
}

void MitchellBestCandidate::compute_min_distances_(std::vector<angem::Point<3,double>> &candidates,
                                                   std::vector<angem::Point<3,double>> &placed,
                                                   std::vector<double> &distances,
                                                   size_t n_placed) const
{
  std::fill( distances.begin(), distances.end(), std::numeric_limits<double>::max() );

  for (size_t i = 0; i < candidates.size(); ++i)
    for (size_t j = 0; j < n_placed; ++j)
      distances[i] = std::min( distances[i], candidates[i].distance( placed[j] ) );
}

}  // end namespace gprs_data
