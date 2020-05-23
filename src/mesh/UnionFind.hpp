#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <unordered_map>
#include <stdexcept>  // invalid_argument
#include <algorithm>  // for_each
#include <numeric>    // iota

namespace mesh {

namespace utils {

using namespace std;

/* Weighted QuickUnion with balancing of the trees
 * Always put the smaller tree lower
 * lookup worst case scenario O(log(n))
 * union O(log(n))
 * Added path compression:
 * access log*(n) -- in practise less than 5 ~ constant
  */
class UnionFind
{
 public:
  UnionFind(const size_t set_size);
  bool connected (const size_t a, const size_t b);
  void merge(const size_t a, const size_t b);
  size_t n_groups() const {return groups.size();}
  size_t find_max(const size_t item) {return largest.find(root(item))->second;}

 protected:
  // size_t root(size_t entry) const;
  // firts value - root
  // second value - depth
  size_t root(size_t entry); // not const because of path compresison
  vector<size_t> id;
  vector<size_t> sizes;
  unordered_set<size_t> groups;
  unordered_map<size_t, size_t> largest;
  friend class Groups;
};


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

} // end namespace utils

}  // end namespace mesh
