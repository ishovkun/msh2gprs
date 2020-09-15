#pragma once
#include "angem/Hexahedron.hpp"
#include <array>

namespace gprs_data {

class UniformCartesianGrid {
 public:
  UniformCartesianGrid(const angem::Point<3,double> &origin,
                       const angem::Point<3,double> &stepping,
                       const std::array<size_t,3> &dims);
  size_t nx() const noexcept { return _dims[0]; }
  size_t ny() const noexcept { return _dims[1]; }
  size_t nz() const noexcept { return _dims[2]; }
  double hx() const noexcept { return _stepping[0]; }
  double hy() const noexcept { return _stepping[1]; }
  double hz() const noexcept { return _stepping[2]; }
  angem::Point<3,double> stepping() const noexcept {return _stepping;}
  double stepping(const size_t dir) const;
  angem::Point<3,double> origin() const noexcept {return _origin;}
  std::array<size_t,3> dims() const noexcept {return _dims;}
  angem::Hexahedron<double> get_bounding_box() const noexcept;

  size_t n_cells() const noexcept;
  size_t n_vertices() const noexcept;
  size_t cell_index(int i, int j, int k) const;
  std::array<int,3> get_ijk(size_t idx) const;
  size_t vertex_index(int i, int j, int k) const;
  bool is_valid_cell(int i, int j, int k) const noexcept;
  angem::Hexahedron<double> get_voxel(size_t idx) const;
  std::vector<size_t> neighbors(size_t search_cell) const;
  bool in_bounds(const angem::Point<3,double> & p) const noexcept;
  void log_stats() const noexcept;
  size_t find_cell(const angem::Point<3,double> & p) const;
  angem::Point<3,double> vertex(size_t v)  const;
  angem::Point<3,double> vertex(int i, int j, int k)  const;
  bool is_valid_vertex(int i, int j, int k) const noexcept;
  std::vector<size_t> neighbors(size_t cell_idx);

 private:
  void add_neighbor_(int i, int j, int k, std::vector<size_t> & dst) const;

  const angem::Point<3,double> _origin;
  const angem::Point<3,double> _stepping;
  const std::array<size_t,3> _dims;
};

}  // end namespace gprs_data
