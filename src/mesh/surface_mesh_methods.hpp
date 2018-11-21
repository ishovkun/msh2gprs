#pragma once

// #include <cstddef>  // std::size_t
#include <vector>

namespace mesh
{

extern const std::size_t MAX_EDGES;

std::size_t estimate_max_edges();


inline std::size_t hash_value(const std::size_t ind1,
                              const std::size_t ind2)
{
  if (ind1 < ind2)
    return MAX_EDGES*ind1 + ind2;
  else
    return MAX_EDGES*ind2 + ind1;
}


inline std::size_t hash_value (const std::pair<std::size_t, std::size_t> & edge)
{
  return hash_value(edge.first, edge.second);
}


inline
std::pair<std::size_t, std::size_t> invert_hash(const std::size_t hash)
{
  std::pair<std::size_t,std::size_t> pair;
  pair.first = (hash - hash % MAX_EDGES) / MAX_EDGES;
  pair.second = hash % MAX_EDGES;
  return pair;
}

}
