#pragma once

#include <vector>
#include <fstream>
#include <iostream>
#include <stack>
#include <limits>

namespace algorithms {

using std::iostream;
using std::fstream;
using std::vector;

class Graph
{
 public:
  Graph(const size_t nv);
  Graph(fstream & in);
  void add(const size_t v1,  const size_t v2);
  const vector<size_t> & adj(const size_t v) const;
  size_t size() const;
  bool hasCycle() const;
  size_t degree(const size_t v) const { return adj(v).size(); }

  friend std::ostream & operator<<(std::ostream & os, const Graph & gr);

 protected:
  vector<vector<size_t>> _adj;
};


}  // end namespace algorithms
