#include <FlowData.hpp>

#include <algorithm> // std::sort
#include <iostream>  // debug


FlowData::FlowData(const std::size_t max_connections)
    :
    max_connections(max_connections)
{}


// void FlowData::reserve_extra(const std::size_t n_elements,
//                              const std::size_t n_connections)
// {
//   const std::size_t old_size = volumes.size();

//   v_neighbors.resize(old_size + n_elements);
//   volumes.resize(old_size + n_elements);
//   poro.resize(old_size + n_elements);
//   depth.resize(old_size + n_elements);

// }



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
  // auto it = element_connection.find(d);
  // if (it == element_connection.end())
  //   throw std::runtime_error("connection apparently does not exist");
  // const std::vector<std::size_t> neighbors = it->second;

  // auto it = element_connection.find(d);
  if (d > v_neighbors.size())
    throw std::runtime_error("connection apparently does not exist");
  const std::vector<std::size_t> neighbors = v_neighbors[d];

  for (const std::size_t neighbor : neighbors)
  {
    const std::size_t dead_conn = connection_index(d, neighbor);
    if (neighbor != u)
    {
      if (!connection_exists(u, neighbor))
      {
        std::size_t new_conn = insert_connection(u, neighbor);
        trans_ij.push_back(trans_ij[dead_conn]);
      }
      else
      {
        const std::size_t conn = connection_index(u, neighbor);
        trans_ij[conn] += trans_ij[dead_conn];
      }
    }
  }

  delete_element(d);
}


void FlowData::clear_connection(const std::size_t ielement,
                                const std::size_t jelement)
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
  // remove connections from v_neighbors
  if (element >= v_neighbors.size())
    throw std::out_of_range("element does not exist: " + std::to_string(element));

  std::vector<std::size_t> neighbors = v_neighbors[element];
  for (const std::size_t neighbor : neighbors)
    clear_connection(neighbor, element);

  v_neighbors.erase(v_neighbors.begin() + element);

  // shift indices
  // in v_neighbors
  for (auto & neighbors : v_neighbors)
    for (std::size_t i=0; i<neighbors.size(); ++i)
      if (neighbors[i] > element)
        neighbors[i]--;

  // rebuild hash in map_connection
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
      throw std::runtime_error("DEBUG: this should be already deleted");
      map_connection.erase(key);
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


bool FlowData::connection_exists(const std::size_t ielement,
                                 const std::size_t jelement) const
{
  const std::size_t hash = hash_value(ielement, jelement);
  auto it = map_connection.find(hash);
  if (it == map_connection.end())
    return false;
  else
    return true;
}


std::size_t FlowData::insert_connection(const std::size_t ielement,
                                        const std::size_t jelement)
{
  const std::size_t hash = hash_value(ielement, jelement);
  const std::size_t conn = map_connection.size();
  if (map_connection.find(hash) != map_connection.end())
    throw std::runtime_error("connection exists");
  map_connection.insert({hash, conn});

  if (std::max(ielement, jelement) >= v_neighbors.size())
  {
    const std::size_t new_size = 2 * std::max(ielement, jelement);
    std::cout << "old size = " << v_neighbors.size() << std::endl;
    std::cout << "new_size = " << new_size << std::endl << std::flush;
    for (const auto & neighbors : v_neighbors)
    {
      // std::cout << "bla-bla" << std::endl;
      std::cout << neighbors.size() << ": ";
      for (const auto & n : neighbors)
        std::cout << n << "\t";
      std::cout << std::endl;
    }

    v_neighbors.resize(new_size);
    std::cout << "resized" << std::endl << std::flush;
  }

  v_neighbors[ielement].push_back(jelement);
  v_neighbors[jelement].push_back(ielement);

  return conn;
}
