#pragma once
#include "SupportRegionsBase.hpp"
#include "discretization/flow/ConnectionData.hpp"
#include "algorithms/EdgeWeightedGraph.hpp"
#include "algorithms/EdgeWeightedDigraph.hpp"
#include <array>
#include <unordered_map>

namespace multiscale {

class SupportRegionsFVMGraph : public SupportRegionsBase  {
 public:
  SupportRegionsFVMGraph(std::vector<size_t> const &partition,
                         algorithms::EdgeWeightedGraph &&connections);
  virtual ~SupportRegionsFVMGraph() = default;

  std::vector<size_t> const & centers() const { return _centers; }
  std::vector<size_t> const & get(size_t block) const {return _support[block];}
  std::vector<size_t> const & get_boundary(size_t block) const { return _support_bnd[block];}
  std::vector<std::vector<size_t>> const & get_edges(size_t block) const {return _support_edges[block];}

 protected:
  size_t find_center_(std::vector<size_t> const  & region, std::vector<size_t> const & bnd) const;

  std::tuple< algorithms::EdgeWeightedDigraph, std::vector<size_t>, std::vector<size_t>>
  build_subgraph_(std::vector<size_t> const & blocks) const;

  std::vector<std::vector<size_t>> find_block_boundaries_() const;
  std::vector<size_t> neighbor_blocks_(std::vector<size_t> const &_bnd) const;
  size_t generate_mapping_(std::vector<size_t> const &blocks, std::vector<size_t> &mapping) const;
  void build_support_region_(std::vector<size_t> blocks, size_t region);
  void modify_edge_weights_();

  algorithms::EdgeWeightedGraph _cons;
  std::vector<std::vector<size_t>> _block_bds;
  std::vector<std::vector<size_t>> _support;
  std::vector<std::vector<size_t>> _support_bnd;
  std::vector<std::vector<std::vector<size_t>>> _support_edges;
};

}  // end namespace multiscale
