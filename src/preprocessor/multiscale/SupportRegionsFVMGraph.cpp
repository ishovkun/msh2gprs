#include "SupportRegionsFVMGraph.hpp"
#include "algorithms/EdgeWeightedGraph.hpp"
#include "algorithms/DijkstraSP.hpp"
#include "MetisInterface.hpp"
#include <numeric>  // accumulate

namespace multiscale {

using discretization::ConnectionData;
using namespace algorithms;

SupportRegionsFVMGraph::SupportRegionsFVMGraph(std::vector<size_t> const &partition,
                                               algorithms::EdgeWeightedGraph &&connections)
    : _partition(partition)
    , _cons(connections)
    , SupportRegionsBase()
    , _blocks(find_cells_in_blocks_())
{
  // find coarse centers
  auto block_bnd = find_block_boundaries_();
  for (size_t coarse = 0; coarse < _blocks.size(); ++coarse) {
    size_t const center = find_center_(_blocks[coarse], block_bnd[coarse]);
    _centers.push_back(center);
  }

  // find support regions
  _support.resize(_blocks.size());
  for (size_t coarse = 0; coarse < _blocks.size(); ++coarse) {
    std::vector<size_t> neighbors = neighbor_blocks_(block_bnd[coarse]);
    std::cout << "building support region for block " << coarse << std::endl;
    auto g = build_support_region_graph_(neighbors);
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

std::vector<std::vector<size_t>> SupportRegionsFVMGraph::find_cells_in_blocks_() const
{
  size_t const n_coarse = *std::max_element(_partition.begin(), _partition.end()) + 1;

  std::vector<std::vector<std::size_t>> cells_in_block(n_coarse);
  for (size_t i = 0; i < _partition.size(); ++i)
    cells_in_block[_partition[i]].push_back(i);
  return cells_in_block;
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

// algorithms::EdgeWeightedGraph SupportRegionsFVMGraph::build_support_region_graph_(std::vector<size_t> const & blocks) const
std::vector<int> SupportRegionsFVMGraph::build_support_region_graph_(std::vector<size_t> const & blocks) const
{
  size_t source = 0;
  // size_t sink = 1;
  // size_t const offset = 2;  // 2 = source + sink
  size_t const offset = 1;
  size_t nv = std::accumulate(blocks.begin(), blocks.end(), offset,
                              [this](size_t cur, size_t block) {
                                return _blocks[block].size() + cur;});

  // build mapping
  std::cout << "build maping (blocks ";
  for (auto block : blocks)
    std::cout << block << " ";
  std::cout << ")" << std::endl;

  std::vector<size_t> mapping(nv);
  std::unordered_map<size_t, size_t> mapping_inv;
  mapping[source] = source;
  // mapping[sink] = sink;
  size_t cnt = 1;
  for (size_t const block : blocks) {
    for (size_t cell : _blocks[block]) {
      mapping[cnt] = cell;
      mapping_inv[cell] = cnt;
      cnt++;
    }
  }
  assert( cnt == nv );

  double sum_weights = std::accumulate( _cons.edges().begin(), _cons.edges().end(),
                                        0.f, [](double cur , auto const & edge) {
                                          return cur + edge.weight();
                                        });
  double const large_weight = 1e6 * sum_weights;
  double const tiny_weight = 0.f;

  // build graph
  std::cout << "build local graph from full graph" << std::endl;
  algorithms::EdgeWeightedGraph g(nv);
  for (size_t const block : blocks) {
    for (size_t cell : _blocks[block]) {
      for (auto const * e : _cons.adj(cell)) {
        size_t const v = e->either();
        size_t const w = e->other(v);
        if ( mapping_inv.count(v) && mapping_inv.count(w) )
        {
          double weight = e->weight();
          if (_partition[mapping_inv[v]] != _partition[mapping_inv[w]])
            weight = large_weight;
          g.add(UndirectedEdge(mapping_inv[v], mapping_inv[w], weight));
        }

      }
    }
  }


  std::cout << "modify local graph" << std::endl;
  // add additional graph edges for the proper partitioning
  // for (size_t coarse : blocks) {
  for (size_t c1 = 0; c1 < blocks.size(); ++c1) {
    size_t const coarse = blocks[c1];
    size_t const v = _centers[coarse];
    size_t const v_mapped = mapping_inv[v];
    g.add(UndirectedEdge(source, v_mapped, large_weight));
    // for (size_t c2 = c1+1; c2 < blocks.size(); ++c2) {
    //   size_t coarse2 = blocks[c2];
    //   size_t const w = _centers[coarse];
    //   size_t const w_mapped = mapping_inv[w];
    //   g.add(UndirectedEdge(v_mapped, w_mapped, large_weight));
    // }


    // auto cutting_edge = g.adj( v_mapped ).front();
    // cutting_edge->set_weight(large_weight);

    // size_t const w = cutting_edge->other(v_mapped);
    // g.add(UndirectedEdge(w, sink, tiny_weight));
  }

  std::cout << "invoke metis" << std::endl;
  // partition into 2 regions: support and no-support
  auto part = multiscale::MetisInterface::partition(g, 2);
  std::vector<int> values(_partition.size(), 0);
  std::cout << "partition: ";
  for (size_t i = offset; i < part.size(); ++i)
  {
    values[mapping[i]] = part[i] + 1;
    // std::cout << mapping[i] << " " << part[i] + 1 << "\n";
  }
  std::cout << std::endl;

  std::cout << "done" << std::endl;
  // exit(0);

  return values;
  // return g;
}


}  // end namespace multiscale
