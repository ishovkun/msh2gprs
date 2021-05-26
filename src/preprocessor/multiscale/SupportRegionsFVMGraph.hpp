#pragma once
#include "SupportRegionsBase.hpp"
#include "discretization/flow/ConnectionData.hpp"
#include "algorithms/EdgeWeightedGraph.hpp"
#include <vector>
#include <array>

namespace multiscale {

class SupportRegionsFVMGraph : public SupportRegionsBase  {
 public:
  SupportRegionsFVMGraph(std::vector<size_t> const &partition,
                         algorithms::EdgeWeightedGraph &&connections);
  virtual ~SupportRegionsFVMGraph() = default;

  std::vector<size_t> const & centers() const { return _centers; }
  size_t size() const { return _blocks.size();}
  std::vector<int> const & get(size_t block) const {return _support[block];}

 protected:
  std::vector<std::vector<size_t>> find_cells_in_blocks_() const;
  size_t find_center_(std::vector<size_t> const  & region, std::vector<size_t> const & bnd) const;
  algorithms::EdgeWeightedGraph build_graph_(std::vector<discretization::ConnectionData> const &connections) const;
  std::vector<std::vector<size_t>> find_block_boundaries_() const;
  std::vector<size_t> neighbor_blocks_(std::vector<size_t> const &_bnd) const;
  // algorithms::EdgeWeightedGraph build_support_region_graph_(std::vector<size_t> const & blocks) const;
  std::vector<int> build_support_region_graph_(std::vector<size_t> const & blocks) const;

  std::vector<size_t> const & _partition;
  algorithms::EdgeWeightedGraph _cons;
  std::vector<std::vector<size_t>> _blocks;
  std::vector<size_t> _centers;
  std::vector<std::vector<int>> _support;
};

}  // end namespace multiscale
