#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <unordered_map>
#include <stdexcept>  // invalid_argument
#include <algorithm>  // for_each
#include <numeric>    // iota


namespace algorithms {

using namespace std;

/* Weighted QuickUnion with balancing of the trees
 * Always put the smaller tree lower
 * lookup worst case scenario O(log(n))
 * union O(log(n))
 * Added path compression:
 * access log*(n) -- in practise less than 5 ~ linear
  */
class UnionFind
{
 public:
  UnionFind(const size_t set_size);
  bool connected (const size_t a, const size_t b);
  void merge(const size_t a, const size_t b);
  size_t n_groups() const {return groups.size();}
  size_t find_max(const size_t item) {return largest.find(root(item))->second;}
  inline size_t group(const size_t item) {return root(item);}

 protected:
  // size_t root(size_t entry) const;
  // firts value - root
  // second value - depth
  size_t root(size_t entry); // not const because of path compresison
  vector<size_t> id;
  vector<size_t> sizes;
  unordered_set<size_t> groups;
  unordered_map<size_t, size_t> largest;
};

}
