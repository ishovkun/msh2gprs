#include "SupportRegionsLaplaceMod.hpp"
#include "algorithms/EdgeWeightedGraph.hpp"
#include "ShapeFunctionSolver.hpp"

namespace multiscale {

using namespace algorithms;

void create_or_update_edge(EdgeWeightedGraph & g, size_t u, size_t v, double weight) {
  auto adj = g.adj(u);
  auto it = adj.begin();

  for (; it != adj.end(); ++it)
    if ( (*it)->other( u ) == v)
      break;

  if (it == adj.end()) g.add(UndirectedEdge(u, v, weight));
  else (*it)->set_weight((*it)->weight() + weight);
}

EdgeWeightedGraph build_block_graph(size_t nb, std::vector<size_t> const & partition,
                                    std::vector<discretization::ConnectionData> const & cons)
{
  EdgeWeightedGraph g(nb);
  for (auto const & con : cons) {
    size_t const u = partition[con.elements[0]];
    size_t const v = partition[con.elements[1]];
    if ( u != v) {
      create_or_update_edge(g, u, v, std::fabs(con.coefficients[0]));
    }
  }
  return g;
}

SupportRegionsLaplaceMod::SupportRegionsLaplaceMod(std::vector<size_t> const &partition, gprs_data::SimData &data)
    : _data(data)
    , SupportRegionsBase(partition)
{
  for (size_t coarse = 0; coarse < _blocks.size(); ++coarse) {
    size_t const center = find_center_(_blocks[coarse]);
    _centers.push_back(center);
  }

  auto const block_cons = build_block_graph(_blocks.size(), _partition, _data.flow_connection_data);

  _support.resize(_blocks.size());
  for (size_t u = 0; u < _blocks.size(); ++u)
  {
    std::cout << "building region for block " << u << std::endl;
    std::vector<size_t> bnd;
    std::vector<double> bnd_values;
    for (auto * e : block_cons.adj(u)) {
      bnd.push_back(_centers[e->other(u)]);
      bnd_values.push_back(0.5);
    }
    size_t source = _centers[u];
    bnd.push_back(source);
    bnd_values.push_back(1.f);

    std::vector<size_t> region;
    for (auto u : _blocks[u])
      region.push_back(u);
    for (auto * e : block_cons.adj(u))
      for (auto v : _blocks[e->other(u)])
        region.push_back(v);

    ShapeFunctionSolver solver( source, region, bnd, bnd_values, _data );

    _support[u].resize( _partition.size(), 0.f);
    auto soln = solver.solution();
    for (size_t i = 0; i < region.size(); ++i)
      _support[u][region[i]] = soln[i];
  }
}

size_t SupportRegionsLaplaceMod::find_center_(std::vector<size_t> const  & region) const
{
  angem::Point<3,double> c;
  double volume = 0;
  for (auto v : region) {
    c += _data.cv_data[v].center * _data.cv_data[v].volume;
    volume += _data.cv_data[v].volume;
  }

  c /= volume;
  size_t ans = region.front();
  for (auto v : region) {
    if ( _data.cv_data[v].center.distance(c) < _data.cv_data[ans].center.distance(c) )
      ans = v;
  }
  return ans;
}

}  // end namespace multiscale
