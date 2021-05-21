#include "SupportRegionsFVMGraph.hpp"
#include "EdgeWeightedDigraph.hpp"
#include "DijkstraSP.hpp"

namespace multiscale {

using discretization::ConnectionData;
using hash_algorithms::ConnectionMap;
using namespace algorithms;

SupportRegionsFVMGraph::SupportRegionsFVMGraph(std::vector<size_t> const &partition,
                                               algorithms::EdgeWeightedGraph &&connections)
    : _partition(partition), _cons(connections), SupportRegionsBase(),
      _blocks(find_cells_in_blocks_())
    // , _block_connections(build_block_connections_())
{
  auto block_bnd = find_block_boundaries_();
  for (size_t coarse = 0; coarse < _blocks.size(); ++coarse) {
    size_t const center = find_center_(_blocks[coarse], block_bnd[coarse]);
    std::cout << "center = " << center << std::endl;
    _centers.push_back(center);
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
      std::cout << "searching from " << bnd[i] << std::endl;
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

ConnectionMap<std::vector<std::array<size_t,2>>> SupportRegionsFVMGraph::build_block_connections_() const
{
  ConnectionMap<std::vector<std::array<size_t,2>>> block_connection;
  for (auto const & edge: _cons.edges())
  {
    size_t const v = edge.either();
    size_t const w = edge.other(v);
    size_t const c1 = _partition[v];
    size_t const c2 = _partition[w];

    if (c1 != c2)
    {
      if ( !block_connection.contains(c1, c2) )
        block_connection.insert(c1, c2);
      block_connection.get_data(c1, c2).push_back({v, w});
    }
  }
  return block_connection;
}

}  // end namespace multiscale
