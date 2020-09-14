#include "GridIntersectionSearcher.hpp"

namespace gprs_data {

GridIntersectionSearcher::GridIntersectionSearcher(const mesh::Mesh & grid)
    : _grid(grid), _mapper(grid)
{}


}  // end namespace gprs_data
