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

std::tuple< algorithms::EdgeWeightedDigraph, std::vector<size_t>, std::vector<size_t>>
SupportRegionsFVMGraph::build_subgraph_(std::vector<size_t> const & blocks) const
{
  std::vector<size_t> mapping_inv;
  size_t const nv = generate_mapping_(blocks, mapping_inv);
  EdgeWeightedDigraph g(nv);
  for (size_t const block : blocks) {
    add_edges(g, _blocks[block], _cons, mapping_inv);
  }

  std::vector<size_t> mapping(nv);
  size_t v = 0;
  for (size_t const b : blocks)
    for (size_t const cell : _blocks[b])
      mapping[v++] = cell;

  std::vector<size_t> centers;
  for (size_t const b : blocks)
    centers.push_back( mapping_inv[_centers[b]] );

  return std::tie( g, mapping, centers );
}

size_t SupportRegionsFVMGraph::find_center_(std::vector<size_t> const  & region, std::vector<size_t> const & bnd) const
{
  // size_t const n = region.size();
  std::fill(_subgraph_mask.begin(), _subgraph_mask.end(), false);
  for (size_t const v : region)
    _subgraph_mask[v] = true;

  // // map between cell indices and local indices within block
  // std::unordered_map<size_t,size_t> mapping;
  // for (size_t iv = 0; iv < n; ++iv) {
  //   mapping[region[iv]] = iv;
  // }

  // graph that consists of only cells in the current block
  // EdgeWeightedDigraph g( n );
  // for (size_t iv = 0; iv < n; ++iv) {
  //   size_t const v = region[iv];
  //   assert( iv == mapping[v] );
  //   for (auto const * const e: _cons.adj(v)) {
  //     size_t const w = e->other(v);
  //     if ( _partition[v] == _partition[w] ) {  // from the same coarse block
  //       g.add(DirectedEdge(iv, mapping[w], e->weight()));
  //       g.add(DirectedEdge(mapping[w], iv, e->weight()));
  //     }
  //   }
  // }

  // find block center as the vertex with the longest path from the boundary
  std::vector<double> farthest(region.size(), std::numeric_limits<double>::max());
  static constexpr double fraction = 0.5;
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
      // DijkstraSP path(g, mapping[bnd[i]]);
      // for (size_t j = 0; j < n; ++j)
      //   farthest[j] = std::min(path.distanceTo(mapping[ region[j] ]), farthest[j]);
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

size_t SupportRegionsFVMGraph::generate_mapping_(std::vector<size_t> const &blocks,
                                                 std::vector<size_t> &mapping) const
{
  size_t const nv = std::accumulate(blocks.begin(), blocks.end(), 0,
                                    [this](size_t cur, size_t block) {
                                      return _blocks[block].size() + cur;});
  mapping.assign(_cons.n_vertices(), nv);
  size_t v = 0;
  for (size_t const block : blocks)
    for (size_t cell : _blocks[block])
      mapping[cell] = v++;
  return nv;
}

std::vector<double> compute_distances(DijkstraSP const &path,
                                      std::vector<size_t> const & centers)
{
  // find distances of neighbor block centers from the current block center
  std::vector<double> dist(centers.size(), 0);
  for (size_t ib = 0; ib < centers.size(); ++ib)
    dist[ib] = path.distanceTo(centers[ib]);
  return dist;
}

std::vector<angem::Point<3,double>> get_triangle(double d12, double d13, double d23)
{
  double const cosa = (d13*d13 + d12*d12 - d23*d23) / (2.f * d13 * d12);
  double const sina = std::sqrt(1.f - cosa*cosa);
  return {
    { 0.f        , 0.f, 0.f },
    { d12        , 0.f, 0.f },
    { d13 * cosa , d13*sina, 0.f}
  };
}

std::vector<bool> flag_vertices(std::vector<size_t> const & mapping,
                                std::vector<DijkstraSP> const &paths,
                                std::vector<double> const & dist,
                                std::vector<size_t> const &partition,
                                std::vector<size_t> const &centers)
{
  std:vector<bool> flags(mapping.size(), false);
  for (size_t v = 0; v < mapping.size(); ++v) {
    size_t const block = partition[v];

    if (block == paths.size() - 1) {
      flags[v] = true;
    }
    else {
      double threshold = 0;
      double norm = 0;

      std::vector<double> dist(paths.size()-1);
      for (size_t ib = 0; ib < paths.size() - 1; ++ib)
        dist[ib] = paths[ib].distanceTo(v);
      std::vector<size_t> idx(dist.size());
      std::iota(idx.begin(), idx.end(), 0);
      std::sort(idx.begin(), idx.end(), [&dist](size_t i1, size_t i2) {return dist[i1] < dist[i2];});

      // size_t b1 = idx[0], b2 = idx[1];

      // double w1 = paths.back().distanceTo(centers[b1]);
      // double w2 = paths.back().distanceTo(centers[b2]);
      // get_triangle(double d12, double d13, double d23);

      // closest
      // threshold = paths.back().distanceTo(centers[b1]);

      // linear
      // if (block != b1 && block != b2)
      //   b2 = block;
      // double val1 = paths.back().distanceTo(centers[b1]);
      // double val2 = paths.back().distanceTo(centers[b2]);
      // double d1 = dist[b1];
      // double d2 = dist[b2];
      // double xi = d1 / (d1 + d2);
      // threshold = val1*(1.f-xi) + val2*xi;

      // linear geometric
      // double b1b2 = paths[b1].distanceTo(centers[b2]);
      // double b1v = dist[b1];
      // double b2v = dist[b2];
      // double cosa = (b1v*b1v + b1b2*b1b2 - b2v*b2v) / (2.f * b1v * b1b2);
      // double b1p = b1v * cosa;
      // double xi = b1p / b1b2;
      // if (xi < 0)      threshold = w1;
      // else if (xi > 1) threshold = w2; // i doubt it will come to that since b1 is closest
      // else             threshold = w1 * (1.f - xi) + w2*xi;

      // FEM triangle
      // double b1b2 = paths[b1].distanceTo(centers[b2]);
      // std::vector<angem::Point<3,double>> vertices(3);
      // vertices[0] = { 0, 0, 0 };
      // vertices[1] = { b1b2, 0, 0 };
      // vertices[2] = { b1b2, 0, 0 };

      // IDW 3 points
      for (size_t i = 0; i < std::min(idx.size() - 1, (size_t)3); ++i) {
        size_t ib = idx[i];
        double value = paths.back().distanceTo(centers[ib]);
        double d = paths[ib].distanceTo(v);
        double weight = 1.f / std::pow(d, 2.0);
        threshold += weight * value;
        norm += weight;

        if (d == 0) {
          threshold = value;
          norm = 1.f;
          break;
        }
      }
      threshold /= norm;

      if (paths.back().distanceTo(v) < threshold)
        flags[v] = true;
    }
  }

  return flags;
}

std::vector<size_t> find_boundary(algorithms::EdgeWeightedDigraph const & g,
                                  std::vector<bool> const &flags,
                                  std::vector<size_t> const & mapping,
                                  EdgeWeightedGraph const & global)
{
  std::vector<size_t> bnd;
  for (auto const & e : g.edges())
    if ((int)flags[e.to()] + (int)flags[e.from()] == 1) {
      if (flags[e.to()]) bnd.push_back(mapping[ e.to() ]);
      else bnd.push_back(mapping[ e.from() ]);
    }

  for (size_t v = 0; v < g.n_vertices(); ++v)
    if (flags[v] && g.adj(v).size() < global.adj(mapping[v]).size())
      bnd.push_back(mapping[v]);

  return bnd;
}

void SupportRegionsFVMGraph::build_support_region_(std::vector<size_t> blocks, size_t region)
{
  blocks.push_back(region);
  // build digraph in order to find shortest paths in it
  auto const [g, mapping, sources] = build_subgraph_(blocks);
  // build shortest paths starting from block centers
  std::vector<DijkstraSP> paths;
  for (size_t ib = 0; ib < blocks.size(); ++ib)
    paths.emplace_back(g, sources[ib]);
  // find distances of neighbor block centers from the current block center
  std::vector<double> const dist = compute_distances(paths.back(), sources);
  // convert global partition vector to include only local blocks
  std::vector<size_t> local(_blocks.size(), blocks.size());
  for (size_t ib = 0; ib < blocks.size(); ++ib)
    local[blocks[ib]] = ib;

  // we don't need target region in the blocks vector any more
  blocks.pop_back();

  std::vector<size_t> part(g.n_vertices());
  for (size_t v = 0; v < g.n_vertices(); ++v)
    part[v] = local[ _partition[mapping[v]] ];

  // flag vertices whether they are in suppport or not
  std::vector<bool> flags = flag_vertices(mapping, paths, dist, part, sources);

  _support_bnd[region] = find_boundary(g, flags, mapping, _cons);
  for (size_t v = 0; v < flags.size(); ++v)
    if (flags[v]) _support[region].push_back(mapping[v]);

  _support_edges.resize(_blocks.size());
  for (size_t ib = 0; ib < blocks.size(); ++ib) {
    _support_edges[region].emplace_back();
    auto &support_edge = _support_edges[region].back();
    for (auto * edge : paths.back().pathTo( sources[ib] )) {
      support_edge.push_back(mapping[edge->to()]);
    }
  }
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
