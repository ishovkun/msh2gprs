#pragma once
#include "SupportRegionsBase.hpp"
#include "discretization/flow/ConnectionData.hpp"
#include "ConnectionMap.hpp"
#include "EdgeWeightedGraph.hpp"
#include <vector>
#include <array>

namespace multiscale {

class SupportRegionsFVMGraph : public SupportRegionsBase  {
 public:
  SupportRegionsFVMGraph(std::vector<size_t> const &partition,
                         algorithms::EdgeWeightedGraph &&connections);
  virtual ~SupportRegionsFVMGraph() = default;

  std::vector<size_t> const & centers() const { return _centers; }

 protected:
  std::vector<std::vector<size_t>> find_cells_in_blocks_() const;
  hash_algorithms::ConnectionMap<std::vector<std::array<size_t,2>>> build_block_connections_() const;
  size_t find_center_(std::vector<size_t> const  & region, std::vector<size_t> const & bnd) const;
  algorithms::EdgeWeightedGraph build_graph_(std::vector<discretization::ConnectionData> const &connections) const;
  std::vector<std::vector<size_t>> find_block_boundaries_() const;


  std::vector<size_t> const & _partition;
  algorithms::EdgeWeightedGraph _cons;
  std::vector<std::vector<size_t>> _blocks;
  hash_algorithms::ConnectionMap<std::vector<std::array<size_t,2>>> _block_connections;
  std::vector<size_t> _centers;
};

}  // end namespace multiscale
