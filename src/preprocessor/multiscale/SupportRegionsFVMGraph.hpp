#pragma once
#include "SupportRegionsBase.hpp"
#include "discretization/flow/ConnectionData.hpp"
#include "algorithms/EdgeWeightedGraph.hpp"
#include <array>
#include <unordered_map>

namespace multiscale {

class SupportRegionsFVMGraph : public SupportRegionsBase  {
 public:
  SupportRegionsFVMGraph(std::vector<size_t> const &partition,
                         algorithms::EdgeWeightedGraph &&connections);
  virtual ~SupportRegionsFVMGraph() = default;

  std::vector<size_t> const & centers() const { return _centers; }
  // std::vector<int> const & get(size_t block) const {return _support[block];}
  std::vector<double> const & get(size_t block) const {return _support[block];}

 protected:
  size_t find_center_(std::vector<size_t> const  & region, std::vector<size_t> const & bnd) const;
  algorithms::EdgeWeightedGraph build_graph_(std::vector<discretization::ConnectionData> const &connections) const;
  std::vector<std::vector<size_t>> find_block_boundaries_() const;
  std::vector<size_t> neighbor_blocks_(std::vector<size_t> const &_bnd) const;
  // algorithms::EdgeWeightedGraph build_support_region_graph_(std::vector<size_t> const & blocks) const;
  std::unordered_map<size_t, size_t> generate_mapping_(std::vector<size_t> const &blocks) const;
  std::vector<double> build_support_region_graph2_(std::vector<size_t> blocks, size_t region) const;

  algorithms::EdgeWeightedGraph _cons;
  // std::vector<std::vector<int>> _support;
  std::vector<std::vector<double>> _support;
  std::vector<std::vector<size_t>> _block_bnd;
};

}  // end namespace multiscale
