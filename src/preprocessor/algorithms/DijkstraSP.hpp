#pragma once
#include "SingleSourceShortestPath.hpp"
#include "IndexMinPQ.hpp"

namespace algorithms {

/* Uses E log(V) to find the path:
 * we process an order of log(V) edges, and use
 * a priority queue with decrease key capability that stores vertices,
 * so each operation there takes log V. */
class DijkstraSP : public SingleSourceShortestPath {
 public:
  DijkstraSP(const EdgeWeightedDigraph & gr, const size_t source)
      : SingleSourceShortestPath(gr, source),
        _pq(_gr.n_vertices())
  {
    _edge_to.resize(_gr.n_vertices());
    _dist_to.resize(_gr.n_vertices(), std::numeric_limits<double>::max());
    _dist_to[_source] = 0;
    _pq.enqueue( _source, 0.0 );

    while (!_pq.empty())
    {
      const size_t v = _pq.dequeue();
      for (const auto & edge : _gr.adj(v))
        relax_(edge);
    }
  }

  virtual ~DijkstraSP() = default;

 protected:
  virtual void relax_(const DirectedEdge & edge) override
  {
    const double new_dist = _dist_to[edge.from()] + edge.weight();
    if ( new_dist < _dist_to[edge.to()] )
    {
      _dist_to[edge.to()] = new_dist;
      _edge_to[edge.to()] = edge;
      // priority is gonna get higher since new_dist < _dist_to[edge.to()]
      if (_pq.contains(edge.to())) _pq.setPriority(edge.to(), edge.weight());
      else                         _pq.enqueue(edge.to(), edge.weight());
    }
  }

 private:
  IndexMinPQ _pq;
};

}  // end namespace algorithms
