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
      for (const auto * edge : _gr.adj(v))
        relax_(*edge);
    }
  }

  virtual ~DijkstraSP() = default;

 protected:
  virtual void relax_(const DirectedEdge & edge) override
  {
    size_t const u = edge.from(), v = edge.to();
    const double new_dist = _dist_to[u] + edge.weight();
    if ( new_dist < _dist_to[v] )
    {
      _dist_to[v] = new_dist;
      _edge_to[v] = edge;
      // priority increases since new distance is shorter
      if (_pq.contains(edge.to())) _pq.setPriority(v, _dist_to[v]);
      else                         _pq.enqueue(v, _dist_to[v]);
    }
  }

 private:
  IndexMinPQ _pq;
};

}  // end namespace algorithms
