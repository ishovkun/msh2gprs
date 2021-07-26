#include "SupportRegionsFVMGraph.hpp"
#include "algorithms/EdgeWeightedGraph.hpp"
#include "algorithms/DijkstraSP.hpp"
#include "algorithms/DijkstraSubgraph.hpp"
#include "MetisInterface.hpp"
#include <numeric>  // accumulate

namespace multiscale {

using discretization::ConnectionData;
using namespace algorithms;

void homogenize(algorithms::EdgeWeightedGraph & g) {
  for (auto & e : g.edges())
    e.set_weight(1.f);
}

SupportRegionsFVMGraph::SupportRegionsFVMGraph(std::vector<size_t> const &partition,
                                               algorithms::EdgeWeightedGraph &&connections)
    : _cons(connections)
    , SupportRegionsBase(partition)
    , _subgraph_mask(_cons.n_vertices(), false)
    , _block_bds(find_block_boundaries_())
{
  modify_edge_weights_();

  // find coarse centers
  for (size_t coarse = 0; coarse < _blocks.size(); ++coarse) {
    size_t const center = find_center_(_blocks[coarse], _block_bds[coarse]);
    _centers.push_back(center);
  }

  // homogenize(_cons);

  // find support regions
  _support_bnd.resize(_blocks.size());
  _support.resize(_blocks.size());
  for (size_t coarse = 0; coarse < _blocks.size(); ++coarse) {
    std::vector<size_t> neighbors = neighbor_blocks_(_block_bds[coarse]);
    std::cout << "building support region for block " << coarse << std::endl;
    build_support_region_(neighbors, coarse);
  }
}

void add_edges(EdgeWeightedDigraph & g,
               std::vector<size_t> const & cells,
               EdgeWeightedGraph const & global,
               std::vector<size_t> const & mapping)
{
  for (size_t u : cells) {
    for (auto const * e : global.adj(u)) {
      size_t const v = e->other(u);
      if ( mapping[u] < g.n_vertices() && mapping[v] < g.n_vertices() )
        g.add(DirectedEdge(mapping[u], mapping[v], e->weight()));
    }
  }
}

size_t SupportRegionsFVMGraph::find_center_(std::vector<size_t> const  & region, std::vector<size_t> const & bnd) const
{
  std::fill(_subgraph_mask.begin(), _subgraph_mask.end(), false);
  for (size_t const v : region)
    _subgraph_mask[v] = true;

  // find block center as the vertex with the longest path from the boundary
  std::vector<double> farthest(region.size(), std::numeric_limits<double>::max());
  static constexpr double fraction = 0.9;
  static constexpr size_t min_selected = 2;
  size_t n_selected = std::max(min_selected, (size_t)(bnd.size() * fraction)); // actual number of items extracted
  std::vector<bool> selected(bnd.size(), false);
  std::fill(selected.begin(), selected.begin() + n_selected, true);
  std::random_shuffle( selected.begin(), selected.end() );
  for (size_t i = 0; i < bnd.size(); ++i)
    if ( selected[i] ) {
      DijkstraSubgraph path(bnd[i], _cons, _subgraph_mask);
      for (size_t j = 0; j <  region.size(); ++j)
        farthest[j] = std::min(path.distanceTo(region[j]), farthest[j]);
    }

  size_t const center = std::distance(farthest.begin(), std::max_element(farthest.begin(), farthest.end()));
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

std::vector<size_t> find_boundary(algorithms::EdgeWeightedDigraph const & g,
                                  std::vector<bool> const &flags,
                                  std::vector<size_t> const & mapping,
                                  EdgeWeightedGraph const & global)
{
  // std::vector<size_t> bnd;
  std::unordered_set<size_t> bnd;
  for (auto const & e : g.edges())
    if ((int)flags[e.to()] + (int)flags[e.from()] == 1) {
      if (flags[e.to()]) bnd.insert(mapping[ e.to() ]);
      // else bnd.push_back(mapping[ e.from() ]);
    }

  for (size_t v = 0; v < g.n_vertices(); ++v)
    if (flags[v] && g.adj(v).size() < global.adj(mapping[v]).size())
      bnd.insert(mapping[v]);

  return std::vector<size_t>(bnd.cbegin(), bnd.cend());
}

bool in_support(size_t u,
                std::vector<size_t> const & blocks,
                std::vector<DijkstraSubgraph> const & paths,
                std::vector<size_t> const & centers,
                std::vector<double> & dist,  // avoid creating all the time
                std::vector<size_t> & idx)   // avoid creating all the time
{
  // sort blocks by distance from vertex u
  for (size_t source = 0; source < blocks.size(); ++source)
    dist[source] = paths[source].distanceTo(u);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&dist](size_t i1, size_t i2) {return dist[i1] < dist[i2];});

  // inverse distance weighting (idw)
  static constexpr size_t n_idw_blocks = 3;  // pick only several closest blocks
  double threshold = 0.f;
  double norm = 0.f;
  for (size_t i = 0; i < std::min(idx.size(), n_idw_blocks); ++i) {
    size_t const source = idx[i];
    double const value = paths.back().distanceTo(centers[blocks[source]]); // distance from region to neighbor center
    double const d = paths[source].distanceTo(u);
    double weight = 1.f / std::pow(d, 2.0);
    threshold += weight * value;
    norm += weight;

    if (d == 0) {  // if distance is zero then we are in the block center
      threshold = value;
      norm = 1.f;
      break;
    }
  }

  threshold /= norm;
  return paths.back().distanceTo(u) < threshold;
}

void SupportRegionsFVMGraph::build_support_region_(std::vector<size_t> blocks, size_t region)
{
  std::fill( _subgraph_mask.begin(), _subgraph_mask.end(), false );
  blocks.push_back(region);
  for (size_t block : blocks)
    for (size_t v : _blocks[block])
      _subgraph_mask[v] = true;

  std::vector<DijkstraSubgraph> paths;
  for (size_t const block : blocks)
    paths.emplace_back( _centers[block], _cons, _subgraph_mask );

  // we don't need target region in the blocks vector any more
  blocks.pop_back();

  // full support region
  // cells that belong to region are always in support
  std::unordered_set<size_t> support;
  for (size_t const v : _blocks[region])
    support.insert(v);
  std::vector<double> dist(blocks.size());
  std::vector<size_t> idx(blocks.size());
  for (size_t ib = 0; ib < blocks.size(); ++ib) {
    size_t const b = blocks[ib];
    for (size_t const u : _blocks[b]) {
      if (in_support(u, blocks, paths, _centers, dist, idx))
        support.insert(u);
    }
  }
  // build boundary and clean it from support
  std::unordered_set<size_t> bnd;
  for (size_t const v : support)
    for (auto const * e : _cons.adj(v))
      if ( !support.count(e->other(v)) )
        bnd.insert( v );
  for (size_t const v : bnd)
    support.erase( support.find(v) );

  _support[region] = std::move(std::vector<size_t>(support.begin(), support.end()));
  _support_bnd[region] = std::move( std::vector<size_t>(bnd.begin(), bnd.end()) );
}

void SupportRegionsFVMGraph::modify_edge_weights_()
{
  double min_nonzero_weight = std::numeric_limits<double>::max();
  for (auto const & e : _cons.edges())
    if (e.weight() != 0.f)
      min_nonzero_weight = std::min(min_nonzero_weight, e.weight());

  double max_weight = 1.f / min_nonzero_weight + 1.f;
  for (auto & e : _cons.edges())
    if (e.weight() == 0.f)
      e.set_weight(max_weight);
    else e.set_weight( 1.f / e.weight() );
}

}  // end namespace multiscale
