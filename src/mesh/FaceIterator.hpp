#pragma once

#include <uint256/uint256_t.h>
#include <angem/Point.hpp>
#include <vector>
#include <unordered_map>

namespace mesh
{
using Point = angem::Point<3,double>;

class FaceIterator
{
 public:
  FaceIterator(const std::size_t iface,
               std::unordered_map<uint256_t, std::vector<std::size_t>> & map_faces,
               std::unordered_map<uint256_t, int> map_physical_faces);
  Point center() const;
  Point normal() const;

 private:
  std::size_t iface;
  std::unordered_map<uint256_t, std::vector<std::size_t>> & map_faces;
  std::unordered_map<uint256_t, int> & map_physical_faces;
};

}
