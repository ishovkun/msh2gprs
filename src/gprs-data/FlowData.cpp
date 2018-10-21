#include <FlowData.hpp>

#include <algorithm> // std::sort
#include <iostream>  // debug


FlowData::FlowData(const std::size_t max_connections)
    :
    max_connections(max_connections)
{}


void FlowData::merge_elements(const std::size_t updated_element,
                              const std::size_t merged_element)
{
  const std::size_t
      u = updated_element,
      d = merged_element;

  // update cell data
  const double v0 = volumes[u];
  const double v1 = volumes[d];
  volumes[u] += v1;
  poro[u] += (v0*poro[u] + v1*poro[d]) / (v0 + v1);
  depth[u] += (v0*depth[u] + v1*depth[d]) / (v0 + v1);

  for (std::size_t i=0; i<custom_data[u].size(); ++i)
    custom_data[u][i] = (v0*custom_data[u][i] + v1*custom_data[d][i]) / (v0 + v1);

  // update face data
  // connections and transes with other elements
  auto it = element_connection.find(d);
  if (it == element_connection.end())
    throw std::runtime_error("connection apparently does not exist");
  const std::vector<std::size_t> neighbors = it->second;

  for (const std::size_t neighbor : neighbors)
  {
    std::cout << "checking dead" << std::endl;
    const std::size_t dead_conn = connection_index(d, neighbor);
    std::cout << "dead done" << std::endl;
    if (neighbor != u)
    {
      std::size_t new_conn = insert_connection(u, neighbor);
      trans_ij.push_back(trans_ij[dead_conn]);
    }
    clear_connection(d, neighbor);
  }

  delete_element(d);
}


void FlowData::clear_connection(const std::size_t ielement,
                                const std::size_t jelement)
{
  auto it = element_connection.find(ielement);
  std::size_t counter = 0;
  for (auto e : it->second)
  {
    if (e == jelement)
      it->second.erase(it->second.begin() + counter);
    counter++;
  }
  if (it->second.size() == 0)
    element_connection.erase(it);

  it = element_connection.find(jelement);
  counter = 0;
  for (auto e : it->second)
  {
    if (e == ielement)
      it->second.erase(it->second.begin() + counter);
    counter++;
  }
  if (it->second.size() == 0)
    element_connection.erase(it);

  // clear other container
  const std::size_t hash = hash_value(ielement, jelement);
  auto it_dead_conn = map_connection.find(hash);
  std::size_t dead_conn = it_dead_conn->second;
  map_connection.erase(it_dead_conn);
  trans_ij.erase(trans_ij.begin() + dead_conn);

  // shift connection indices
  for (auto iter = map_connection.begin();
       iter != map_connection.end(); ++iter)
    if (iter->second > dead_conn)
      iter->second--;
}


void FlowData::delete_element(const std::size_t element)
{
  // rebuild keys (expensive)
  std::vector<std::size_t> keys;
  keys.reserve(map_connection.size());
  for (auto it : map_connection)
    keys.push_back(it.first);
  std::sort(keys.begin(), keys.end());

  for (const auto & key : keys)
  {
    std::size_t conn = map_connection[key];
    std::pair<std::size_t,std::size_t> elements = invert_hash(key);
    if (elements.first > element and elements.second < element)
    {
      std::size_t new_hash = hash_value(elements.first-1, elements.second);
      map_connection.erase(key);
      map_connection.insert({ new_hash, conn });
    }
    else if (elements.first < element and elements.second > element)
    {
      std::size_t new_hash = hash_value(elements.first, elements.second-1);
      map_connection.erase(key);
      map_connection.insert({ new_hash, conn });
    }
    else if (elements.first > element and elements.second > element)
    {
      std::size_t new_hash = hash_value(elements.first-1, elements.second-1);
      map_connection.erase(key);
      map_connection.insert({ new_hash, conn });
    }
    else if (elements.first == element or elements.second == element)
    {
      map_connection.erase(key);
    }
  }

  // update neighbor map
  keys.clear();
  keys.reserve(element_connection.size());
  for (auto it : element_connection)
    keys.push_back(it.first);
  std::sort(keys.begin(), keys.end());

  for (const auto & key : keys)
  {
    if (key == element)
    {
      // const auto neighbors = element_connection[key];
      element_connection.erase(key);
    }
    else
    {
      for (auto & ielement : element_connection[key])
      {
        if (ielement > element)
          ielement--;
      }
    }
  }

  // delete cell data
  volumes.erase(volumes.begin() + element);
  poro.erase(poro.begin() + element);
  depth.erase(depth.begin() + element);
  custom_data.erase(custom_data.begin() + element);

}


std::size_t FlowData::connection_index(const std::size_t ielement,
                                       const std::size_t jelement) const
{
  const std::size_t hash = hash_value(ielement, jelement);
  auto it = map_connection.find(hash);
  if (it == map_connection.end())
    throw std::runtime_error("connection does not exist");
  return it->second;
}
