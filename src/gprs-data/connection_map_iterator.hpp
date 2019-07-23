#pragma once

#include <vector>
#include <unordered_map>

namespace hash_algorithms
{

template <typename DataType>
class connection_map_iterator : public std::iterator<std::forward_iterator_tag, DataType>
{
 public:
  // constructor
  connection_map_iterator(std::unordered_map<std::size_t, std::size_t>::iterator it,
                          std::vector<DataType>                                & data,
                          const std::size_t                                      max_elements);
  // ITERATOR OPERATORS:
  // comparison operator
  bool operator==(const connection_map_iterator<DataType> & other);
  // negative comparision operator
  bool operator!=(const connection_map_iterator<DataType> & other);
  connection_map_iterator operator++() {map_it++; return *this;}
  // data access operator
  DataType & operator*() const {return data[map_it->second];}
  // data access operator
  DataType * operator->() {return &data[map_it->second];}
  // get element indices from connection iterator
  inline std::pair<std::size_t,std::size_t> elements() const
  {return invert_hash(map_it->first);}
  std::size_t connection_index() const {return map_it->second;}

 protected:
  std::pair<std::size_t,std::size_t> invert_hash(const std::size_t hash_value) const;

 private:
  const std::size_t max_elements;
  std::unordered_map<std::size_t, std::size_t>::iterator map_it;
  std::vector<DataType> & data;
};


template <typename DataType>
connection_map_iterator<DataType>::
connection_map_iterator(std::unordered_map<std::size_t, std::size_t>::iterator it,
                        std::vector<DataType>                                & data,
                        const std::size_t                                      max_elements)
    : map_it(it), data(data), max_elements(max_elements)
{}


template <typename DataType>
inline bool connection_map_iterator<DataType>::
operator==(const connection_map_iterator<DataType> & other)
{
  return (map_it == other.map_it);
}


template <typename DataType>
inline bool connection_map_iterator<DataType>::
operator!=(const connection_map_iterator<DataType> & other)
{
  return (map_it != other.map_it);
}


template <typename DataType>
inline std::pair<std::size_t,std::size_t> connection_map_iterator<DataType>::
invert_hash(const std::size_t hash) const
{
  std::pair<std::size_t,std::size_t> pair;
  pair.first = (hash - hash % max_elements) / max_elements;
  pair.second = hash % max_elements;
  return pair;
}

/*  ------------------------------------- Const-iterator ----------------- */
template <typename DataType>
class connection_map_const_iterator : public std::iterator<std::forward_iterator_tag, DataType>
{
 public:
  // constructor
  connection_map_const_iterator(std::unordered_map<std::size_t, std::size_t>::const_iterator it,
                                const std::vector<DataType>                          & data,
                                const std::size_t                                      max_elements);
  // ITERATOR OPERATORS:
  // comparison operator
  bool operator==(const connection_map_const_iterator<DataType> & other) const;
  // negative comparision operator
  bool operator!=(const connection_map_const_iterator<DataType> & other) const;
  connection_map_const_iterator operator++() {map_it++; return *this;}
  // data access operator
  DataType & operator*() const {return data[map_it->second];}
  // get element indices from connection iterator
  inline std::pair<std::size_t,std::size_t> elements() const
  {return invert_hash(map_it->first);}

 protected:
  std::pair<std::size_t,std::size_t> invert_hash(const std::size_t hash_value) const;

 private:
  const std::size_t max_elements;
  std::unordered_map<std::size_t, std::size_t>::const_iterator map_it;
  const std::vector<DataType> & data;
};


template <typename DataType>
connection_map_const_iterator<DataType>::
connection_map_const_iterator(std::unordered_map<std::size_t, std::size_t>::const_iterator it,
                              const std::vector<DataType>                          & data,
                              const std::size_t                                      max_elements)
    : map_it(it), data(data), max_elements(max_elements)
{}


template <typename DataType>
inline bool connection_map_const_iterator<DataType>::
operator==(const connection_map_const_iterator<DataType> & other) const
{
  return (map_it == other.map_it);
}


template <typename DataType>
inline bool connection_map_const_iterator<DataType>::
operator!=(const connection_map_const_iterator<DataType> & other) const
{
  return (map_it != other.map_it);
}


template <typename DataType>
inline std::pair<std::size_t,std::size_t> connection_map_const_iterator<DataType>::
invert_hash(const std::size_t hash) const
{
  std::pair<std::size_t,std::size_t> pair;
  pair.first = (hash - hash % max_elements) / max_elements;
  pair.second = hash % max_elements;
  return pair;
}


}
