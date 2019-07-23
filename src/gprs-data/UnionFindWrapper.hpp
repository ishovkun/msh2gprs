#pragma once

#include "UnionFind.hpp"
#include <unordered_map>
#include <memory>  // std::shared_prt

namespace algorithms
{

template <typename T>
class UnionFindWrapper
{
 public:
  UnionFindWrapper();
  void insert(const T & item);
  void finalize();
  void merge(const T & item1, const T & item2);
  const std::unordered_map<T, size_t> & items() const {return storage;}
  size_t group(const T & item) const;

 private:
  std::unordered_map<T, size_t> storage;
  std::unique_ptr<UnionFind> p_uf;
  bool is_finalized;
};


template <typename T>
UnionFindWrapper<T>::UnionFindWrapper()
    :
    is_finalized(false)
{}


template <typename T>
void UnionFindWrapper<T>::finalize()
{
  if (is_finalized)
    throw std::runtime_error("union_find already finalized");

  p_uf = std::make_unique<UnionFind>(storage.size());
  is_finalized = true;
}


template <typename T>
void UnionFindWrapper<T>::insert(const T & item)
{
  if (is_finalized)
    throw std::runtime_error("union_find already finalized");

  storage.insert({item, storage.size()});
}


template <typename T>
void UnionFindWrapper<T>::merge(const T & item1,
                                const T & item2)
{
  assert(is_finalized);
  p_uf->merge(storage[item1], storage[item2]);
}


template <typename T>
size_t UnionFindWrapper<T>::group(const T & item) const
{
  if (!is_finalized)
    throw std::runtime_error("union_find not finalized");
  const auto it = storage.find(item);
  if (it == storage.end()) throw std::invalid_argument("item does not exist");
  return p_uf->group(it->second);
}


}  // end namespace
