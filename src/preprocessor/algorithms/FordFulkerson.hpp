#pragma once
#include "FlowNetwork.hpp"
#include <limits>
#include <queue>

namespace algorithms {

class FordFulkerson {
 public:
  FordFulkerson(FlowNetwork & gr, size_t source, size_t sink)
      : _gr(gr), _source(source), _sink(sink), _marked(gr.n_vertices()), _edge_to(gr.n_vertices())
  {
    _total_flow = 0;
    _niter = 0;
    while (hasAugmentingPath())
    {
      // find minimum bottleneck
      double bottle = std::numeric_limits<double>::max();
      for (size_t v = _sink; v != _source; v = _edge_to[v]->other(v))
        bottle = std::min(bottle, _edge_to[v]->residualCapacityTo(v));
      // augment the path
      for (size_t v = _sink; v != _source; v = _edge_to[v]->other(v))
        _edge_to[v]->addResidualFlowTo(v, bottle);

      _total_flow += bottle;
      _niter++;
    }
  }

  // get total value of flow
  double total_flow() const noexcept {return _total_flow;}
  // is reachable from source in residual network
  bool inCut(size_t v) const {return _marked[v];}

  // number of augmenting iterations
  size_t niter() const noexcept {return _niter;}

  virtual ~FordFulkerson() = default;

 private:
  bool hasAugmentingPath()
  {
    std::fill(_marked.begin(), _marked.end(), false);
    std::fill(_edge_to.begin(), _edge_to.end(), nullptr);
    // BFS to find the path
    _marked[_source] = true;
    std::queue<size_t> q;
    q.push(_source);
    while (!q.empty())
    {
      const size_t v = q.front();
      q.pop();
      for (const auto * edge : _gr.adj(v))
      {
        const size_t w = edge->other(v);
        if (!_marked[w] && edge->residualCapacityTo(w) > 0)
        {
          _edge_to[w] = const_cast<FlowEdge*>(edge);
          _marked[w] = true;
          q.push(w);
        }
      }
    }
    return _marked[_sink];
  }


  FlowNetwork & _gr;
  const size_t _source, _sink;
  std::vector<bool> _marked;
  std::vector<FlowEdge*> _edge_to;
  double _total_flow = 0;
  size_t _niter = 0;
};


}  // end namespace algorithms
