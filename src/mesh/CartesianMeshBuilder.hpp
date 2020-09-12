#pragma once

#include "Mesh.hpp"

namespace mesh {

struct CartesianMeshParameters
{
  std::vector<double> dx = {1.0};
  std::vector<double> dy = {1.0};
  std::vector<double> dz = {1.0};
  angem::Point<3, double> origin;
};

class CartesianMeshBuilder {
 public:
  CartesianMeshBuilder(const CartesianMeshParameters & data);
  size_t nx() const noexcept;
  size_t ny() const noexcept;
  size_t nz() const noexcept;
  size_t nvx() const noexcept;
  size_t nvy() const noexcept;
  size_t nvz() const noexcept;
  size_t cell_index(size_t i, size_t j, size_t k) const noexcept;
  size_t vertex_index(size_t i, size_t j, size_t k) const noexcept;
  size_t n_cells() const noexcept;
  size_t n_vertices() const noexcept;
  operator Mesh() const;

 private:
  void setup_vertices_(Mesh & grid) const;
  void setup_cells_(Mesh & grid) const;

  CartesianMeshParameters _data;
};

}  // end namespace mesh
