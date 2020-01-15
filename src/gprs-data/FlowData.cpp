#include <FlowData.hpp>

#include <algorithm> // std::sort
#include <iostream>  // debug


namespace flow
{

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

  const double v0 = cells[u].volume;
  const double v1 = cells[d].volume;
  cells[u].volume += v1;
  cells[u].porosity = (v0*cells[u].porosity + v1*cells[d].porosity) / (v0 + v1);
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
        flow::FaceData new_face = insert_connection(d, neighbor);
        flow::FaceData modified_face = get_connection(u, neighbor);
        new_face.transmissibility = modified_face.transmissibility;
        new_face.thermal_conductivity = modified_face.thermal_conductivity;
        new_face.conType = modified_face.conType;
        std::size_t new_face_conN = modified_face.conCV.size();
        new_face.conCV.resize(new_face_conN);
        new_face.conTr.resize(new_face_conN);
        new_face.conArea.resize(new_face_conN);
        new_face.conPerm.resize(new_face_conN);
        for(std::size_t m=0; m < new_face_conN; m++){
            new_face.zVolumeFactor.resize(new_face_conN);
            if(modified_face.conCV[m] == d){
                new_face.conCV[m] = u;
            } else{
                new_face.conCV[m] = modified_face.conCV[m];
            }
            new_face.conTr[m] = modified_face.conTr[m];
            new_face.conArea[m] = modified_face.conArea[m];
            new_face.conPerm[m] = modified_face.conPerm[m];
            new_face.zVolumeFactor[m] = modified_face.zVolumeFactor[m];
        }
      }
      else
      {
        flow::FaceData modified_face = get_connection(u, neighbor);
        std::size_t modified_face_conN = modified_face.conCV.size();
        modified_face.transmissibility += dead_connection.transmissibility;
        modified_face.thermal_conductivity += dead_connection.thermal_conductivity;
        for(std::size_t m=0; m < modified_face_conN; m++){
            modified_face.conTr[m] += dead_connection.conTr[m];
            modified_face.conArea[m] += dead_connection.conArea[m];
            modified_face.zVolumeFactor[m] = (v0*modified_face.zVolumeFactor[m] + v1*dead_connection.zVolumeFactor[m]) / (v0 + v1);
            modified_face.conPerm[m] = modified_face.conTr[m]*modified_face.zVolumeFactor[m]/modified_face.conArea[m];
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
