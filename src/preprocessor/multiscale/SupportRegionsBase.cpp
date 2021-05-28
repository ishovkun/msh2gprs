#include "SupportRegionsBase.hpp"
#include <algorithm>

namespace multiscale {

SupportRegionsBase::SupportRegionsBase(std::vector<size_t> const &partition)
    :_partition(partition), _blocks(find_cells_in_blocks_())
{}

std::vector<std::vector<size_t>> SupportRegionsBase::find_cells_in_blocks_() const
{
  size_t const n_coarse = *std::max_element(_partition.begin(), _partition.end()) + 1;

  std::vector<std::vector<std::size_t>> cells_in_block(n_coarse);
  for (size_t i = 0; i < _partition.size(); ++i)
    cells_in_block[_partition[i]].push_back(i);
  return cells_in_block;
}


}  // end namespace multiscale
