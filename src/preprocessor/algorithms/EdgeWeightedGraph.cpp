#include "EdgeWeightedGraph.hpp"

namespace algorithms {

using namespace std;

EdgeWeightedGraph::EdgeWeightedGraph(const size_t nv)
    : _adj(nv)
{}

void EdgeWeightedGraph::add(const UndirectedEdge & e)
{
  const size_t ie = _edges.size();
  _edges.push_back(e);
  const size_t v = e.either();
  const size_t w = e.other(v);
  _adj[v].push_back(ie);
  _adj[w].push_back(ie);
}

std::ostream & operator<<(ostream & os, const EdgeWeightedGraph & gr)
{
  os << gr.n_vertices() << "\n";
  os << gr.n_edges() << "\n";
  for (const auto & e : gr.edges())
    os << e.either() << " " << e.other(e.either()) << " " << e.weight() << "\n";
  return os;
}

std::vector<UndirectedEdge const *> EdgeWeightedGraph::adj(const size_t v) const
{
  std::vector<UndirectedEdge const*> result;
  for (const size_t ie : _adj[v])
  {
    UndirectedEdge const & e = _edges[ie];
    result.push_back(&e);
  }
  return result;
}

std::vector<UndirectedEdge*> EdgeWeightedGraph::adj(const size_t v)
{
  std::vector<UndirectedEdge*> result;
  for (const size_t ie : _adj[v])
  {
    UndirectedEdge & e = _edges[ie];
    result.push_back(&e);
  }
  return result;
}

bool EdgeWeightedGraph::has_edge(size_t u, size_t v) const
{
  for (const size_t ie : _adj[u]) {
    UndirectedEdge const & e = _edges[ie];
    if ( e.other(u) == v )
      return true;
  }

  return false;
}


}  // end namespace algorithms
