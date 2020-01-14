#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>  // debug


namespace flow
{


struct CellData
{
  double volume;
  double porosity;
  double depth;
  std::vector<double> custom;
};


struct FaceData
{
  double transmissibility; // geometric part of transmissibility of each connection
  double thermal_conductivity; // thermal conductivity of each connection

  std::size_t conType; // connection type (1:M-M, 2:M-F, 3:F-F)
  std::vector<double> conCV; // Index of control volumes (elements) involved in each face
  // Geometric part of a transmissibility of each control volume
  // = (abs perm) * (face area) / (distance between CV center and a face)
  std::vector<double> conTr;
  std::vector<double> conArea; // Face area [m2]
  std::vector<double> zVolumeFactor; // aperture [m]
  std::vector<double> conPerm; // Absolute permeability of each control volume [md]
};


class FlowData
{
 public:
  FlowData(const std::size_t max_connections = 1e10);
  void reserve_extra(const std::size_t n_elements,
                     const std::size_t n_connections);
  // returns the connection index
  FaceData & insert_connection(const std::size_t ielement,
                             const std::size_t jelement);
  // throws std::out_of_range if connection does not exist
  FaceData & get_connection(const std::size_t ielement,
                            const std::size_t jelement);
  bool connection_exists(const std::size_t ielement,
                         const std::size_t jelement) const;
  // get two elements from hash value
  std::pair<std::size_t,std::size_t>
  invert_hash(const std::size_t hash) const;
  void merge_elements(const std::size_t updated_element,
                      const std::size_t merged_element);
  void clear_connection(const std::size_t ielement,
                        const std::size_t jelement);
  void delete_element(const std::size_t element);


 private:
  std::size_t hash_value(const std::size_t ielement,
                         const std::size_t jelement) const;

 public:
  std::unordered_map<std::size_t, FaceData> map_connection;
  std::vector<std::vector<std::size_t>> v_neighbors;

  std::vector<CellData> cells;
  std::vector<FaceData> faces;
  std::vector<std::string>         custom_names;

 private:
  std::size_t max_connections;
};


// inline


inline
std::size_t FlowData::hash_value(const std::size_t ind1,
                                 const std::size_t ind2) const
{
  if (ind1 < ind2)
    return max_connections * ind1 + ind2;
  else
    return max_connections * ind2 + ind1;
}


inline
std::pair<std::size_t,std::size_t>
FlowData::invert_hash(const std::size_t hash) const
{
  std::pair<std::size_t,std::size_t> pair;
  pair.first = (hash - hash % max_connections) / max_connections;
  pair.second = hash % max_connections;

  if (map_connection.find(hash) == map_connection.end())
    throw std::runtime_error("element does not exist");

  return pair;
}

}
