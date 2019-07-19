#include "UnionFind.hpp"

namespace algorithms {

UnionFind::UnionFind(const size_t set_size)
    :id(set_size), sizes(set_size)
{
  std::iota(id.begin(), id.end(), 0);
  std::fill(sizes.begin(), sizes.end(), 1);
  for (size_t i=0; i < set_size; i++)
  {
    groups.insert(i);
    largest.insert({i, i});
  }
}


size_t UnionFind::root(size_t entry) //const
{
  // without path compression
  //  while (entry != id[entry]) entry = id[entry];
  //  with path compression
  while (entry != id[entry])
  {
    id[entry] = id[id[entry]];
    entry = id[entry];
  }
  return entry;
}


bool UnionFind::connected(const size_t a, const size_t b)
{
  return root(a) == root(b);
}


void UnionFind::merge(const size_t a, const size_t b)
{
  const auto root_a = root(a);
  const auto root_b = root(b);
  if (root_a == root_b) return;
  if (sizes[root_a] >= sizes[root_b])
  {
    id[root_b] = root_a;
    sizes[root_a] += sizes[root_b];
    groups.erase(groups.find(root_b));
    const size_t max = std::max(largest[root_a], largest[root_b]);
    largest.erase(largest.find(root_b));
    largest[root_a] = max;
  }
  else
  {
    id[root_a] = root_b;
    sizes[root_b] += sizes[root_a];
    groups.erase(groups.find(root_a));
    const size_t max = std::max(largest[root_a], largest[root_b]);
    largest.erase(largest.find(root_a));
    largest[root_b] = max;
  }
}


}
