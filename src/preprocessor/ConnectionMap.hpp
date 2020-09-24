#pragma once

#include <connection_map_iterator.hpp>
#include <vector>
#include <unordered_map>
#include <algorithm> // std::sort


namespace hash_algorithms
{

// empty struct that allows connectionmap not to
// waste space on connection data
struct empty {};

template <typename DataType>
class ConnectionMap
{
 public:
  ConnectionMap(const std::size_t max_elements = 1e10);

  connection_map_iterator<DataType> begin() {
    return connection_map_iterator<DataType>(connections.begin(), m_data, max_elements);
  }
  connection_map_iterator<DataType> end() {
    return connection_map_iterator<DataType>(connections.end(), m_data, max_elements);
  }
  connection_map_const_iterator<DataType> begin() const {
    return connection_map_const_iterator<DataType>(connections.begin(), m_data, max_elements);
  }
  connection_map_const_iterator<DataType> end() const {
    return connection_map_const_iterator<DataType>(connections.end(), m_data, max_elements);
  }

  // get number of connections
  std::size_t size() const {return connections.size();}

  inline DataType & get_data(const std::size_t connection_index)
  {
    assert( connection_index < m_data.size() );
    return m_data[connection_index];
  }
  inline const DataType & get_data(const std::size_t connection_index) const {return m_data[connection_index];}
  DataType & get_data(const std::size_t ielement, const std::size_t jelement);
  const DataType & get_data(const std::size_t ielement, const std::size_t jelement) const;

  // creates a new connection and returns the connection index
  std::size_t insert(const std::size_t ielement, const std::size_t jelement);
  // remove all the data about the element in internal storage
  void delete_element(const std::size_t ielement);
  // throws std::out_of_range if connection does not exist
  std::size_t index(const std::size_t ielement, const std::size_t jelement) const;
  // returns true if the connection between elements exists
  bool contains(const std::size_t ielement, const std::size_t jelement) const;
  // returns true if the connection between elements exists
  bool contains(const std::pair<size_t,size_t> & elements) const;
  // get neighbors of the connection map element
  const std::vector<std::size_t> & get_neighbors(std::size_t ielement) const;
  // delete a connection between elements from the map
  // does not clear the storage
  void remove(const std::size_t ielement, const std::size_t jelement);

 private:
  // delete a connection between elements from the map
  void clear_(const std::size_t ielement, const std::size_t jelement);

  std::size_t hash_value(const std::size_t ind1,
                         const std::size_t ind2) const;

  std::pair<std::size_t,std::size_t> invert_hash(const std::size_t hash) const;
  void merge_elements(const std::size_t updated_element,
                      const std::size_t merged_element);

  const std::size_t max_elements;
  std::unordered_map<std::size_t, std::size_t> connections;
  std::vector<std::vector<std::size_t>> v_neighbors;

  std::vector<DataType> m_data;
};


template <typename DataType>
ConnectionMap<DataType>::ConnectionMap(const std::size_t max_elements)
    : max_elements(max_elements)
{}


template <typename DataType>
inline std::size_t
ConnectionMap<DataType>::hash_value(const std::size_t ind1,
                                    const std::size_t ind2) const
{
  if (ind1 < ind2)
    return max_elements * ind1 + ind2;
  else
    return max_elements * ind2 + ind1;
}


template <typename DataType>
inline std::pair<std::size_t,std::size_t>
ConnectionMap<DataType>::invert_hash(const std::size_t hash) const
{
  std::pair<std::size_t,std::size_t> pair;
  pair.first = (hash - hash % max_elements) / max_elements;
  pair.second = hash % max_elements;

  if (connections.find(hash) == connections.end())
    throw std::runtime_error("element does not exist");

  return pair;
}

template <typename DataType>
void ConnectionMap<DataType>::remove(const std::size_t ielement, const std::size_t jelement)
{
  const std::size_t hash = hash_value(ielement, jelement);
  auto it = connections.find(hash);
  if (it != connections.end())
    connections.erase(it);
  else throw std::invalid_argument("connection " + std::to_string(ielement) + " " +
                                   std::to_string(jelement) + " does not exist");
}

template <typename DataType>
void ConnectionMap<DataType>::clear_(const std::size_t ielement, const std::size_t jelement)
{
  // update neighbors vector
  if (ielement >= v_neighbors.size())
    throw std::out_of_range("element does not exist: " + std::to_string(ielement));
  if (jelement >= v_neighbors.size())
    throw std::out_of_range("element does not exist: " + std::to_string(jelement));

  for (std::size_t i=0; i<v_neighbors[ielement].size(); ++i)
  {
    if (v_neighbors[ielement][i] == jelement)
    {
      v_neighbors[ielement].erase(v_neighbors[ielement].begin() + i);
      break;
    }
  }
  for (std::size_t i=0; i<v_neighbors[jelement].size(); ++i)
  {
    if (v_neighbors[jelement][i] == ielement)
    {
      v_neighbors[jelement].erase(v_neighbors[jelement].begin() + i);
      break;
    }
  }

  // clear other container
  const std::size_t hash = hash_value(ielement, jelement);
  auto it_dead_conn = connections.find(hash);
  std::size_t dead_conn = it_dead_conn->second;
  connections.erase(it_dead_conn);
  m_data.erase(m_data.begin() + dead_conn);

  // shift connection indices
  for (auto iter = connections.begin();
       iter != connections.end(); ++iter)
    if (iter->second > dead_conn)
      iter->second--;
}


template <typename DataType>
void ConnectionMap<DataType>::delete_element(const std::size_t element)
{
  // remove connections from v_neighbors
  if (element >= v_neighbors.size())
    throw std::out_of_range("element does not exist: " + std::to_string(element));

  std::vector<std::size_t> neighbors = v_neighbors[element];
  for (const std::size_t neighbor : neighbors)
    clear_(neighbor, element);

  v_neighbors.erase(v_neighbors.begin() + element);

  // shift indices
  // in v_neighbors
  for (auto & neighbors : v_neighbors)
    for (std::size_t i=0; i<neighbors.size(); ++i)
      if (neighbors[i] > element)
        neighbors[i]--;

  // rebuild hash in connections
  std::vector<std::size_t> keys;
  keys.reserve(connections.size());
  for (auto it : connections)
    keys.push_back(it.first);
  std::sort(keys.begin(), keys.end());
  for (const auto & key : keys)
  {
    std::size_t conn = connections[key];
    std::pair<std::size_t,std::size_t> elements = invert_hash(key);
    if (elements.first > element and elements.second < element)
    {
      std::size_t new_hash = hash_value(elements.first-1, elements.second);
      connections.erase(key);
      connections.insert({ new_hash, conn });
    }
    else if (elements.first < element and elements.second > element)
    {
      std::size_t new_hash = hash_value(elements.first, elements.second-1);
      connections.erase(key);
      connections.insert({ new_hash, conn });
    }
    else if (elements.first > element and elements.second > element)
    {
      std::size_t new_hash = hash_value(elements.first-1, elements.second-1);
      connections.erase(key);
      connections.insert({ new_hash, conn });
    }
    else if (elements.first == element or elements.second == element)
    {
      throw std::runtime_error("DEBUG: this should be already deleted");
      connections.erase(key);
    }
  }
}


template <typename DataType>
std::size_t ConnectionMap<DataType>::index(const std::size_t ielement,
                                           const std::size_t jelement) const
{
  const std::size_t hash = hash_value(ielement, jelement);
  auto it = connections.find(hash);
  if (it == connections.end())
    throw std::runtime_error("connection " + std::to_string(ielement) + "-"
                             + std::to_string(jelement) + " does not exist");
  return it->second;
}


template <typename DataType>
bool ConnectionMap<DataType>::contains(const std::size_t ielement,
                                       const std::size_t jelement) const
{
  const std::size_t hash = hash_value(ielement, jelement);
  auto it = connections.find(hash);
  if (it == connections.end())
    return false;
  else
    return true;
}

template <typename DataType>
bool ConnectionMap<DataType>::contains(const std::pair<size_t,size_t> & elements) const
{
  return contains(elements.first, elements.second);
}

template <typename DataType>
std::size_t ConnectionMap<DataType>::insert(const std::size_t ielement,
                                            const std::size_t jelement)
{
  const std::size_t hash = hash_value(ielement, jelement);
  const std::size_t conn = connections.size();
  if (connections.find(hash) != connections.end())
    throw std::runtime_error("connection exists");
  connections.insert({hash, conn});

  if (std::max(ielement, jelement) >= v_neighbors.size())
  {
    const std::size_t new_size = 2 * std::max(ielement, jelement);
    v_neighbors.resize(new_size);
  }

  v_neighbors[ielement].push_back(jelement);
  v_neighbors[jelement].push_back(ielement);

  if (sizeof(DataType) > sizeof(empty))  // don't waste space if hash_algorithms::empty
    m_data.resize(m_data.size() + 1);

  return conn;
}


template <typename DataType>
const std::vector<std::size_t> &
ConnectionMap<DataType>::get_neighbors(std::size_t ielement) const
{
  assert(ielement < v_neighbors.size());
  return v_neighbors[ielement];
}


template <typename DataType>
DataType &
ConnectionMap<DataType>::get_data(const std::size_t ielement, const std::size_t jelement)
{
  const size_t con_id = index(ielement, jelement);
  assert( con_id <  m_data.size());
  return m_data[con_id];
}


template <typename DataType>
const DataType &
ConnectionMap<DataType>::get_data(const std::size_t ielement, const std::size_t jelement) const
{
  const size_t con_id = index(ielement, jelement);
  assert( con_id <  m_data.size());
  return m_data[con_id];
}

}  // end namespace

using PureConnectionMap = hash_algorithms::ConnectionMap<hash_algorithms::empty>;
