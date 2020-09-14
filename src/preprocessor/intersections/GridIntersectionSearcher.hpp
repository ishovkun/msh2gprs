#pragma once
#include "mesh/Mesh.hpp"
#include "Mapper.hpp"

namespace gprs_data {

class GridIntersectionSearcher {
 public:
  GridIntersectionSearcher(const mesh::Mesh & grid);
  // bool remove(size_t cell_index);
  // void insert(size_t cell_index);

 private:
  const mesh::Mesh & _grid;
  Mapper _mapper;
  // std::unique_ptr<GeometrySearchTreeNode> _root = nullptr;
};

}  // end namespace gprs_data
