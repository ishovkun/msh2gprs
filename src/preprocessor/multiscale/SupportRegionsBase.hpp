#pragma once
#include <vector>

namespace multiscale {

class SupportRegionsBase {
 public:
  SupportRegionsBase(std::vector<size_t> const &partition);
  virtual ~SupportRegionsBase() = default;

  virtual std::vector<size_t> const & centers() const { return _centers; }
  virtual size_t size() {return _blocks.size();}

 protected:
  std::vector<std::vector<size_t>> find_cells_in_blocks_() const;

  std::vector<size_t> const _partition;
  std::vector<std::vector<size_t>> _blocks;
  std::vector<size_t> _centers;
};

}  // end namespace multiscale
