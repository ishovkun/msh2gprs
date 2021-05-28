#include "SupportRegionsFVMGraph.hpp"
#include "algorithms/EdgeWeightedGraph.hpp"
#include "algorithms/DijkstraSP.hpp"
#include "algorithms/FlowNetwork.hpp"
#include "algorithms/FordFulkerson.hpp"
#include "MetisInterface.hpp"
#include <numeric>  // accumulate

namespace multiscale {

using discretization::ConnectionData;
using namespace algorithms;

SupportRegionsFVMGraph::SupportRegionsFVMGraph(std::vector<size_t> const &partition,
                                               algorithms::EdgeWeightedGraph &&connections)
    : _cons(connections)
    , SupportRegionsBase(partition)
{
  // find coarse centers
  auto block_bnd = find_block_boundaries_();
  for (size_t coarse = 0; coarse < _blocks.size(); ++coarse) {
    size_t const center = find_center_(_blocks[coarse], block_bnd[coarse]);
    _centers.push_back(center);
  }

  // find support regions
  _block_bnd = block_bnd;
  _support.resize(_blocks.size());
  for (size_t coarse = 0; coarse < _blocks.size(); ++coarse) {
    std::vector<size_t> neighbors = neighbor_blocks_(block_bnd[coarse]);
    std::cout << "building support region for block " << coarse << std::endl;
    // neighbors.push_back(coarse);
    // auto g = build_support_region_graph_(neighbors, coarse);
    auto g = build_support_region_graph2_(neighbors, coarse);
    _support[coarse] = g;
  }
}

size_t SupportRegionsFVMGraph::find_center_(std::vector<size_t> const  & region, std::vector<size_t> const & bnd) const
{
  size_t const n = region.size();

  // map between cell indices and local indices within block
  std::unordered_map<size_t,size_t> mapping;
  for (size_t iv = 0; iv < n; ++iv) {
    mapping[region[iv]] = iv;
  }

  // graph that consists of only cells in the current block
  EdgeWeightedDigraph g( n );
  for (size_t iv = 0; iv < n; ++iv) {
    size_t const v = region[iv];
    assert( iv == mapping[v] );
    for (auto const * const e: _cons.adj(v)) {
      size_t const w = e->other(v);
      if ( _partition[v] == _partition[w] ) {  // from the same coarse block
        g.add(DirectedEdge(iv, mapping[w], e->weight()));
        g.add(DirectedEdge(mapping[w], iv, e->weight()));
      }
    }
  }

  // find block center as the vertex with the longest path from the boundary
  std::vector<double> paths(n, std::numeric_limits<double>::max());
  double fraction = 0.5;
  size_t n_selected = std::max((size_t)2, (size_t)(bnd.size() * fraction)); // actual number of items extracted
  std::vector<bool> selected(bnd.size(), false);
  std::fill(selected.begin(), selected.begin() + n_selected, true);
  std::random_shuffle( selected.begin(), selected.end() );
  for (size_t i = 0; i < bnd.size(); ++i)
    if ( selected[i] ) {
      DijkstraSP path(g, mapping[bnd[i]]);
      for (size_t j = 0; j < n; ++j)
        paths[j] = std::min(path.distanceTo(j), paths[j]);
    }

  size_t const center = std::distance(paths.begin(), std::max_element(paths.begin(), paths.end()));
  return region[center];
}

std::vector<std::vector<size_t>> SupportRegionsFVMGraph::find_block_boundaries_() const
{
  std::vector<std::vector<size_t>> block_bds(_blocks.size());
  for (auto const & edge: _cons.edges())
  {
    size_t const v = edge.either();
    size_t const w = edge.other(v);
    size_t const c1 = _partition[v];
    size_t const c2 = _partition[w];
    if (c1 != c2)
    {
      block_bds[c1].push_back(v);
      block_bds[c2].push_back(w);
    }
  }
  return block_bds;
}

std::vector<size_t> SupportRegionsFVMGraph::neighbor_blocks_(std::vector<size_t> const & bnd) const
{
  // TODO: 1- or 2-level bfs to include edge neighbors
  // std::vector<size_t> ans;
  std::unordered_set<size_t> ans;
  for (size_t const v : bnd) {
    for (auto * e : _cons.adj(v)) {
      size_t const w = e->other(v);
      if ( _partition[v] != _partition[w] ) ans.insert(_partition[w]);
    }
  }
  return std::vector(ans.begin(), ans.end());
}

std::unordered_map<size_t, size_t> SupportRegionsFVMGraph::generate_mapping_(std::vector<size_t> const &blocks) const
{
  size_t offset = 0;
  size_t nv = std::accumulate(blocks.begin(), blocks.end(), offset,
                              [this](size_t cur, size_t block) {
                                return _blocks[block].size() + cur;});
  std::unordered_map<size_t, size_t> mapping_inv;
  size_t cnt = offset;
  for (size_t const block : blocks) {
    for (size_t cell : _blocks[block]) {
      mapping_inv[cell] = cnt;
      cnt++;
    }
  }
  assert( cnt == nv );
  return mapping_inv;
}

std::vector<double> SupportRegionsFVMGraph::build_support_region_graph2_(std::vector<size_t> blocks, size_t region) const
{
  blocks.push_back(region);
  auto mapping = generate_mapping_(blocks);
  size_t const nv = mapping.size();

  EdgeWeightedDigraph g( nv );
  for (size_t const block : blocks) {
    for (size_t v : _blocks[block]) {
      for (auto const * e : _cons.adj(v)) {
        size_t const w = e->other(v);
        if ( mapping.count(v) && mapping.count(w) )
          g.add(DirectedEdge(mapping[v], mapping[w], e->weight()));
      }
    }
  }

  std::vector<DijkstraSP> paths;
  for (size_t ib = 0; ib < blocks.size(); ++ib)
  {
    size_t b = blocks[ib];
    size_t source = mapping[ _centers[b] ];
    paths.emplace_back(g, source);
  }

  // find farthest vertex and assign it as sink
  std::vector<double> dist(blocks.size(), 0);
  std::cout << "distances: " << std::endl;
  for (size_t ib = 0; ib < blocks.size(); ++ib)
  {
    size_t b = blocks[ib];
    size_t u = mapping[ _centers[b] ];
    dist[ib] = paths.back().distanceTo(u);
    std::cout << dist[ib] << " ";
  }
  std::cout << std::endl;

  // find normalization coefficients
  std::vector<double> coefs(blocks.size(), 0);
  for (size_t v = 0; v < g.n_vertices(); ++v)
    for (size_t ib = 0; ib < blocks.size(); ++ib)
      coefs[ib] = std::max(coefs[ib], paths[ib].distanceTo(v));
  std::cout << "coefs" << std::endl;
  for (auto c : coefs)
    std::cout << c << " ";
  std::cout << std::endl;


  std::vector<double> values(_partition.size(), 0);
  for (size_t const block : blocks) {
    for (size_t cell : _blocks[block]) {
      size_t const v = mapping[ cell ];

      if (block == region) {
        values[cell] = 1;
      }
      else {
        double threshold = 0;
        double norm = 0;
        for (size_t ib = 0; ib < blocks.size() - 1; ++ib) {
          threshold += dist[ib] * (1.f - paths[ib].distanceTo(v) / coefs[ib]);
          norm += (1.f - paths[ib].distanceTo(v) / coefs[ib]);
        }
        threshold /= norm;

        // values[cell] = threshold;
        values[cell] = (paths.back().distanceTo(v) >= threshold) ? 2 : 1;
      }
    }
  }

  return values;
}

}  // end namespace multiscale
