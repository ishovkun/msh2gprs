#include "Graph.hpp"

namespace algorithms {

Graph::Graph(const size_t nv)
: _adj(nv)
{}

Graph::Graph(fstream & in)
{
  size_t size, v1, v2;
  in >> size;
  _adj.resize(size);
  in >> size;  // n_conn
  for (std::size_t i=0; i<size; ++i)
  {
    in >> v1;
    in >> v2;
    add(v1, v2);
  }
}

void Graph::add(const size_t v1, const size_t v2)
{
  _adj[v1].push_back(v2);
  _adj[v2].push_back(v1);
}

size_t Graph::size() const
{
  return _adj.size();
}

const vector<size_t> & Graph::adj(const size_t v) const
{
  if (v >= size()) throw std::invalid_argument("vertex " + std::to_string(v) + " does not exist");
  return _adj[v];
}

std::ostream & operator<<(std::ostream & os, const Graph & gr)
{
  os << "Graph (" << gr.size() << ")" << std::endl;
  for (std::size_t i=0; i<gr.size(); ++i)
  {
    os << i << ": ";
    for (const size_t v : gr.adj(i))
      os << v << " ";
    os << std::endl;
  }
  return os;
}

bool Graph::hasCycle() const
{
  if (!size()) return true;

  std::stack<size_t> stack;
  stack.push(0);
  vector<bool> marked(size());
  marked[0] = true;
  vector<size_t> edge(size(), std::numeric_limits<size_t> :: max());

  while (!stack.empty())
  {
    const size_t i = stack.top();
    stack.pop();
    for (const size_t j : adj(i))
      if (!marked[j])
      {
        edge[j] = i;
        marked[j] = true;
        stack.push(j);
      }
      else if (edge[i] != j)
      {
        return true;
      }
  }
  return false;
}

}  // end namespace algorithms
