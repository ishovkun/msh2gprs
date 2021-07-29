#include "GeometricPartition.hpp"

namespace multiscale {

using namespace gprs_data;

GeometricPartition::GeometricPartition(std::array<size_t,3> dims, gprs_data::SimData & data)
    : _data(data)
    , _uniform(build_uniform_grid_(dims))
    , _partition(data.flow.cv.size(), 0)
{
  for (size_t i = 0; i < data.flow.cv.size(); ++i)
  {
    auto const & c = data.flow.cv[i].center;
    _partition[i] = _uniform.find_cell(c);
  }
}

UniformCartesianGrid GeometricPartition::build_uniform_grid_(std::array<size_t,3> dims) const
{
  assert( dims[0] > 0 && dims[1] > 0 && dims[2] > 0 && "Coarse grid must have valid dimensions" );

  const double max_double = std::numeric_limits<double>::max();
  const double min_double = std::numeric_limits<double>::lowest();
  angem::Point<3,double> glob_min, glob_max;
  glob_min[0] = glob_min[1] = glob_min[2] = max_double;
  glob_max[0] = glob_max[1] = glob_max[2] = min_double;

  for (auto const & cv : _data.flow.cv)
  {
    for (size_t i = 0; i < 3; ++i)
    {
      glob_min[i] = std::min(glob_min[i], cv.center[i]);
      glob_max[i] = std::max(glob_max[i], cv.center[i]);
    }
  }

  angem::Point<3,double> stepping;
  for (size_t i = 0; i < 3; ++i)
    stepping[i] = (glob_max[i] - glob_min[i]) / dims[i];
  return UniformCartesianGrid(glob_min, stepping, dims);
}


}  // end namespace multiscale
