#pragma once
#include "UnionFind.hpp"
#include <map>

namespace mesh {

namespace utils {

class Groups {
 public:
  Groups(const std::vector<size_t> & entities)
      : _uf(entities.size())
  {
    for (size_t i = 0; i < entities.size(); ++i)
      _entity_index[entities[i]] = i;
  }

  void merge( const size_t e1, const size_t e2 )
  {
    _uf.merge( _entity_index[e1], _entity_index[e2] );
  }

  std::vector<std::vector<size_t>> get()
  {
    std::map<size_t, size_t> root_indices;
    size_t i = 0;
    for (const size_t root : _uf.groups)
      root_indices[root] = i++;

    std::vector<std::vector<size_t>> groups(_uf.n_groups());
    for (const auto & it : _entity_index)
    {
      const size_t root = _uf.root(it.second);
      groups[root_indices[root]].push_back(it.first);
    }

    return groups;
  }

  virtual ~Groups() = default;

 private:
  std::unordered_map<size_t, size_t> _entity_index;
  UnionFind _uf;
};

}  // end namespace utils


}  // end namespace mesh
