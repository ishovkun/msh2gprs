#include <FlowData.hpp>

#include <algorithm> // std::sort
#include <iostream>  // debug


namespace flow
{

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
  // const double v0 = volumes[u];
  // const double v1 = volumes[d];
  // volumes[u] += v1;
  // poro[u] = (v0*poro[u] + v1*poro[d]) / (v0 + v1);
  // depth[u] = (v0*depth[u] + v1*depth[d]) / (v0 + v1);

  // for (std::size_t i=0; i<custom_data[u].size(); ++i)
  //   custom_data[u][i] = (v0*custom_data[u][i] + v1*custom_data[d][i]) / (v0 + v1);

  const double v0 = cells[u].volume;
  const double v1 = cells[d].volume;
  cells[u].volume += v1;
  cells[u].porosity = (v0*cells[u].porosity + v1*cells[d].porosity) / (v0 + v1); // Jaewoo An
  cells[u].depth = (v0*cells[u].depth + v1*cells[d].depth) / (v0 + v1);

  for (std::size_t i=0; i<cells[u].custom.size(); ++i)
    cells[u].custom[i] = (v0*cells[u].custom[i] + v1*cells[d].custom[i]) / (v0 + v1);

  // update face data
  // connections and transes with other elements
  if (d > v_neighbors.size())
    throw std::runtime_error("connection apparently does not exist");
  const std::vector<std::size_t> neighbors = v_neighbors[d];

  for (const std::size_t neighbor : neighbors)
  {
    // const std::size_t dead_conn = connection_index(d, neighbor);
    flow::FaceData dead_connection = get_connection(d, neighbor);
    if (neighbor != u)
    {
      if (!connection_exists(u, neighbor))
      {
        // std::size_t new_conn = insert_connection(u, neighbor);
        // faces.emplace_back();
        // faces.back().transmissibility = faces[dead_conn].transmissibility;
        // trans_ij.push_back(trans_ij[dead_conn]);
        flow::FaceData new_face = insert_connection(d, neighbor);
        flow::FaceData modified_face = get_connection(u, neighbor);
        new_face.transmissibility = modified_face.transmissibility;
        new_face.thermal_conductivity = modified_face.thermal_conductivity;
        new_face.ConType = modified_face.ConType; // Jaewoo An
        new_face.ConN = modified_face.ConN;
        new_face.ConCV.resize(modified_face.ConN);
        new_face.ConTr.resize(modified_face.ConN);
        new_face.ConArea.resize(modified_face.ConN);
        new_face.ConPerm.resize(modified_face.ConN);
        for(std::size_t m=0; m < modified_face.ConN; m++){
            new_face.ZVolumeFactor.resize(modified_face.ConN);
            if(modified_face.ConCV[m] == d){
                new_face.ConCV[m] = u;
            } else{
                new_face.ConCV[m] = modified_face.ConCV[m];
            }
            new_face.ConTr[m] = modified_face.ConTr[m];
            new_face.ConArea[m] = modified_face.ConArea[m];
            new_face.ConPerm[m] = modified_face.ConPerm[m];
            new_face.ZVolumeFactor[m] = modified_face.ZVolumeFactor[m];
        }
      }
      else
      {
        // const std::size_t conn = connection_index(u, neighbor);
        // trans_ij[conn] += trans_ij[dead_conn];
        flow::FaceData modified_face = get_connection(u, neighbor);
        modified_face.transmissibility += dead_connection.transmissibility;
        modified_face.thermal_conductivity += dead_connection.thermal_conductivity;
        for(std::size_t m=0; m < modified_face.ConN; m++){
            modified_face.ConTr[m] += dead_connection.ConTr[m];
            modified_face.ConArea[m] += dead_connection.ConArea[m];
            modified_face.ZVolumeFactor[m] = (v0*modified_face.ZVolumeFactor[m] + v1*dead_connection.ZVolumeFactor[m]) / (v0 + v1);
            modified_face.ConPerm[m] = modified_face.ConTr[m]*modified_face.ZVolumeFactor[m]/modified_face.ConArea[m];
        }
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
  map_connection.erase(it_dead_conn);

  // shift connection indices
  // for (auto iter = map_connection.begin();
  //      iter != map_connection.end(); ++iter)
  //   if (iter->second > dead_conn)
  //     iter->second--;
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
  // for (auto & neighbors : v_neighbors)
  //   for (std::size_t i=0; i<neighbors.size(); ++i)
  //     if (neighbors[i] > element)
  //       neighbors[i]--;

  // rebuild hash in map_connection
  // std::vector<std::size_t> keys;
  // keys.reserve(map_connection.size());
  // for (auto it : map_connection)
    // keys.push_back(it.first);

  // std::sort(keys.begin(), keys.end());
  // for (const auto & key : keys)
  // {
  //   // std::size_t conn = map_connection[key];
  //   std::pair<std::size_t,std::size_t> elements = invert_hash(key);
  //   if (elements.first > element and elements.second < element)
  //   {
  //     std::size_t new_hash = hash_value(elements.first-1, elements.second);
  //     map_connection.erase(key);
  //     map_connection.insert({ new_hash, conn });
  //   }
  //   else if (elements.first < element and elements.second > element)
  //   {
  //     std::size_t new_hash = hash_value(elements.first, elements.second-1);
  //     map_connection.erase(key);
  //     map_connection.insert({ new_hash, conn });
  //   }
  //   else if (elements.first > element and elements.second > element)
  //   {
  //     std::size_t new_hash = hash_value(elements.first-1, elements.second-1);
  //     map_connection.erase(key);
  //     map_connection.insert({ new_hash, conn });
  //   }
  //   else if (elements.first == element or elements.second == element)
  //   {
  //     throw std::runtime_error("DEBUG: this should be already deleted");
  //     map_connection.erase(key);
  //   }
  // }

  // delete cell data
  cells.erase(cells.begin() + element);
}


FaceData & FlowData::get_connection(const std::size_t ielement,
                                    const std::size_t jelement)
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


FaceData & FlowData::insert_connection(const std::size_t ielement,
                                       const std::size_t jelement)
{
  const std::size_t hash = hash_value(ielement, jelement);

  if (map_connection.find(hash) != map_connection.end())
    throw std::runtime_error("connection exists");

  auto it = map_connection.insert({hash, FaceData()});

  if (std::max(ielement, jelement) >= v_neighbors.size())
  {
    const std::size_t new_size = 2 * std::max(ielement, jelement);
    v_neighbors.resize(new_size);
  }

  v_neighbors[ielement].push_back(jelement);
  v_neighbors[jelement].push_back(ielement);

  // return conn;
  return it.first->second;
}

}
