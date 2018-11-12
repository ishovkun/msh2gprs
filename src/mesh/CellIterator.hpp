#pragma once

#include <uint256/uint256_t.h>
#include <angem/Point.hpp>
#include <vector>
#include <unordered_map>

namespace mesh
{
using Point = angem::Point<3,double>;

class CellIterator
{
 public:
  // constructor
  CellIterator(const std::size_t icell,
               std::unordered_map<uint256_t, std::vector<std::size_t>> & map_faces,
               std::vector<int> & shape_ids,
               std::vector<int> & cell_markers);
  // assignment operator
  CellIterator & operator=(const CellIterator & other);
  // comparison
  bool operator==(const CellIterator & other) const;
  bool operator!=(const CellIterator & other) const;
  // GETTERS
  inline std::size_t cell_index() const {return icell;}
  inline int shape_id() const {return shape_ids[icell];}
  inline int marker() const {return cell_markers[icell];}
  std::vector<std::size_t> neighbor_indices() const;
  std::vector<CellIterator> neighbors() const;
  Point center() const;

  // incrementing
  CellIterator & operator++();
  CellIterator & operator--();

 private:
  // std::vector<std::vector<std::size_t>> & cells;  // vertex indices
  std::size_t icell;
  std::unordered_map<uint256_t, std::vector<std::size_t>> & map_faces;
  std::vector<int> & shape_ids;
  std::vector<int> & cell_markers;
};

}
